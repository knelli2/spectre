// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/Cce/WorldtubeBufferUpdater.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/ComplexModalVector.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/SpinWeighted.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/Cce/ExtractionRadius.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/Version.hpp"
#include "NumericalAlgorithms/SpinWeightedSphericalHarmonics/SwshCoefficients.hpp"
#include "NumericalAlgorithms/SpinWeightedSphericalHarmonics/SwshTags.hpp"
#include "NumericalAlgorithms/SpinWeightedSphericalHarmonics/SwshTransform.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/Math.hpp"
#include "Utilities/Numeric.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace Cce {
namespace detail {
std::pair<size_t, size_t> create_span_for_time_value(
    const double time, const size_t pad, const size_t interpolator_length,
    const size_t lower_bound, const size_t upper_bound,
    const DataVector& time_buffer) {
  ASSERT(
      lower_bound < upper_bound,
      "The supplied `lower_bound` is greater than `upper_bound`, which is not "
      "permitted");
  ASSERT(2 * interpolator_length + pad <= upper_bound,
         "The combined `interpolator_length` and `pad` is too large for the "
         "supplied `upper_bound`.\nupper_bound="
             << upper_bound << "\npad=" << pad
             << "\ninterpolator_length=" << interpolator_length);

  size_t range_start = lower_bound;
  size_t range_end = upper_bound;
  while (range_end - range_start > 1) {
    if (time_buffer[(range_start + range_end) / 2] < time) {
      range_start = (range_start + range_end) / 2;
    } else {
      range_end = (range_start + range_end) / 2;
    }
  }
  // always keep the difference between start and end the same, even when
  // the interpolations starts to get worse
  size_t span_start = lower_bound;
  size_t span_end =
      std::min(interpolator_length * 2 + pad + lower_bound, upper_bound);
  if (range_end + interpolator_length + pad > upper_bound) {
    span_start =
        std::max(upper_bound - (interpolator_length * 2 + pad), lower_bound);
    span_end = upper_bound;
  } else if (range_start + 1 > lower_bound + interpolator_length) {
    span_start = range_start - interpolator_length;
    span_end = range_end + interpolator_length + pad - 1;
  }

  return std::make_pair(span_start, span_end);
}

void set_time_buffer_and_lmax(const gsl::not_null<DataVector*> time_buffer,
                              size_t& l_max, const h5::Dat& data,
                              const bool is_real, const bool is_modal_data) {
  const auto data_table_dimensions = data.get_dimensions();
  const Matrix time_matrix =
      data.get_data_subset(std::vector<size_t>{0}, 0, data_table_dimensions[0]);
  *time_buffer = DataVector{data_table_dimensions[0]};

  for (size_t i = 0; i < data_table_dimensions[0]; ++i) {
    (*time_buffer)[i] = time_matrix(i, 0);
  }

  if (is_modal_data) {
    // If the quantitiy is real, then there's no way for us to do a check just
    // from the number of columns alone. If not, then we expect an even number
    // of total modes (both real and imaginary), so the number of columns should
    // be odd with the addition of the time column
    if (not is_real and data_table_dimensions[1] % 2 != 1) {
      ERROR(
          "Dimensions of subfile "
          << data.subfile_path()
          << " are incorrect. Was expecting an odd number of columns (because "
             "of the time), but got "
          << data_table_dimensions[1] << " instead.");
    }

    // Avoid compiler warning
    const size_t l_plus_one_squared = (data_table_dimensions[1] - 1) / 2;
    l_max =
        static_cast<size_t>(sqrt(static_cast<double>(l_plus_one_squared)) - 1);
  } else {
    // No way to check number of columns for nodal data

    // If number of real values N = (l+1)*(2l+1), then
    // l = ( -3 + sqrt(9 + 8 * (N-1)) ) / 4. But the dimensions[1] includes time
    // and both real and complex values, so we have to accound for those
    const size_t discriminant =
        9 + 8 * ((data_table_dimensions[1] - 1) / 2 - 1);
    l_max =
        static_cast<size_t>((sqrt(static_cast<double>(discriminant)) - 3)) / 4;
  }
}

template <int Spin>
void update_buffer_with_modal_data(
    const gsl::not_null<ComplexModalVector*> buffer_to_update,
    const h5::Dat& read_data, const size_t computation_l_max,
    const size_t l_max, const size_t time_span_start,
    const size_t time_span_end, const bool is_modal_data) {
  constexpr bool is_real = Spin == 0;
  size_t number_of_columns = read_data.get_dimensions()[1];
  const size_t result_l_max = std::min(l_max, computation_l_max);
  if (UNLIKELY(buffer_to_update->size() !=
               square(computation_l_max + 1) *
                   (time_span_end - time_span_start))) {
    ERROR("Incorrect storage size for the data to be loaded in.");
  }
  std::vector<size_t> cols(number_of_columns - 1);
  std::iota(cols.begin(), cols.end(), 1);
  auto data_matrix =
      read_data.get_data_subset<std::vector<std::vector<double>>>(
          cols, time_span_start, time_span_end - time_span_start);
  *buffer_to_update = 0.0;

  const size_t expected_dat_modal_size = square(l_max + 1);
  const size_t expected_dat_nodal_size = (l_max + 1) * (2 * l_max + 1);
  ASSERT(
      cols.size() == is_modal_data ? expected_dat_modal_size
                                   : expected_dat_nodal_size,
      "Incorrect number of columns in Dat file. Expected (excluding time "
      "column) "
          << (is_modal_data ? expected_dat_modal_size : expected_dat_nodal_size)
          << ", but got " << cols.size());
  const size_t computation_modal_size = square(computation_l_max + 1);

  for (size_t time_row = 0; time_row < time_span_end - time_span_start;
       ++time_row) {
    // Handle nodal data specially because we have to transform it
    if (not is_modal_data) {
      // NOLINTNEXTLINE
      const DataVector raw_nodal_data{data_matrix[time_row].data() + 1,
                                      cols.size()};
      // Not const so the view is cleaner. Also no way to avoid a copy if we
      // need to transform it because ComplexDataVector memory layout is
      // different than DataVector
      ComplexDataVector complex_nodal_data =
          std::complex<double>(1.0, 0.0) * raw_nodal_data;
      SpinWeighted<ComplexDataVector, Spin> nodal_data_view;
      nodal_data_view.set_data_ref(make_not_null(&complex_nodal_data));

      SpinWeighted<ComplexModalVector, Spin> single_goldberg_mode_view;
      single_goldberg_mode_view.set_data_ref(
          buffer_to_update->data() + time_row * computation_modal_size,
          computation_modal_size);

      // First transform to coefficients using swsh_transform, and then convert
      // libsharp coefficients into modes
      Spectral::Swsh::libsharp_to_goldberg_modes(
          make_not_null(&single_goldberg_mode_view),
          Spectral::Swsh::swsh_transform(l_max, 1, nodal_data_view), l_max);

      // Modes are already put into the correct buffer so we can continue on
      continue;
    }

    for (int l = 0; l <= static_cast<int>(result_l_max); ++l) {
      for (int m = -l; m <= l; ++m) {
        const size_t computation_goldber_index =
            Spectral::Swsh::goldberg_mode_index(computation_l_max,
                                                static_cast<size_t>(l), m) *
                (time_span_end - time_span_start) +
            time_row;
        // NOLINTBEGIN
        if (is_real) {
          if (m == 0) {
            (*buffer_to_update)[computation_goldber_index] =
                std::complex<double>(
                    data_matrix[time_row][static_cast<size_t>(square(l))], 0.0);
          } else {
            const double factor = (m > 0 or abs(m) % 2 == 0) ? 1.0 : -1.0;
            (*buffer_to_update)[computation_goldber_index] =
                factor * std::complex<double>(
                             data_matrix[time_row][static_cast<size_t>(
                                 square(l) + 2 * abs(m) - 1)],
                             sgn(m) * data_matrix[time_row][static_cast<size_t>(
                                          square(l) + 2 * abs(m))]);
          }
        } else {
          (*buffer_to_update)[computation_goldber_index] = std::complex<double>(
              data_matrix[time_row][2 * Spectral::Swsh::goldberg_mode_index(
                                            l_max, static_cast<size_t>(l), m)],
              data_matrix[time_row][2 * Spectral::Swsh::goldberg_mode_index(
                                            l_max, static_cast<size_t>(l), m) +
                                    1]);
        }
        // NOLINTEND
      }
    }
  }
}

template <typename InputTags>
double update_buffers_for_time(
    const gsl::not_null<Variables<InputTags>*> buffers,
    const gsl::not_null<size_t*> time_span_start,
    const gsl::not_null<size_t*> time_span_end, const double time,
    const size_t computation_l_max, const size_t l_max,
    const size_t interpolator_length, const size_t buffer_depth,
    const DataVector& time_buffer,
    const tuples::tagged_tuple_from_typelist<
        db::wrap_tags_in<Tags::detail::InputDataSet, InputTags>>& dataset_names,
    const h5::H5File<h5::AccessType::ReadOnly>& cce_data_file,
    const bool is_modal_data) {
  if (*time_span_end >= time_buffer.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (*time_span_end > interpolator_length and
      time_buffer[*time_span_end - interpolator_length] > time) {
    // the next time an update will be required
    return time_buffer[*time_span_end - interpolator_length + 1];
  }
  // find the time spans that are needed
  auto new_span_pair = detail::create_span_for_time_value(
      time, buffer_depth, interpolator_length, 0, time_buffer.size(),
      time_buffer);
  *time_span_start = new_span_pair.first;
  *time_span_end = new_span_pair.second;
  // load the desired time spans into the buffers
  tmpl::for_each<InputTags>([&](auto tag_v) {
    using tag = typename decltype(tag_v)::type;
    update_buffer_with_modal_data<tag::type::type::spin>(
        make_not_null(&get(get<tag>(*buffers)).data()),
        cce_data_file.get<h5::Dat>(
            "/" + get<Tags::detail::InputDataSet<tag>>(dataset_names)),
        computation_l_max, l_max, *time_span_start, *time_span_end,
        is_modal_data);
    cce_data_file.close_current_object();
  });
  // the next time an update will be required
  return time_buffer[std::min(*time_span_end - interpolator_length + 1,
                              time_buffer.size() - 1)];
}

}  // namespace detail

MetricWorldtubeH5BufferUpdater::MetricWorldtubeH5BufferUpdater(
    const std::string& cce_data_filename,
    const std::optional<double> extraction_radius, const bool file_is_from_spec,
    const bool is_modal_data)
    : is_modal_data_(is_modal_data),
      cce_data_file_{cce_data_filename},
      filename_{cce_data_filename},
      file_is_from_spec_(file_is_from_spec) {
  set_names<ComplexModalVector>();
  set_names<ComplexDataVector>();

  // 'VersionHist' is a feature written by SpEC to indicate the details of the
  // file format. This line determines whether or not the radial derivatives
  // require renormalization based on whether the SpEC version that produced it
  // was an old one that had a particular normalization bug
  has_version_history_ = cce_data_file_.exists<h5::Version>("/VersionHist");

  extraction_radius_ =
      Cce::get_extraction_radius(cce_data_filename, extraction_radius, true)
          .value();

  detail::set_time_buffer_and_lmax(make_not_null(&time_buffer_), l_max_,
                                   cce_data_file_.get<h5::Dat>("/Lapse"), false,
                                   is_modal_data_);
  cce_data_file_.close_current_object();
}

template <typename T>
void MetricWorldtubeH5BufferUpdater::set_names() {
  get<Tags::detail::InputDataSet<Tags::detail::SpatialMetric<T>>>(
      dataset_names_) = "/g";
  get<Tags::detail::InputDataSet<
      Tags::detail::Dr<Tags::detail::SpatialMetric<T>>>>(dataset_names_) =
      "/Drg";
  get<Tags::detail::InputDataSet<::Tags::dt<Tags::detail::SpatialMetric<T>>>>(
      dataset_names_) = "/Dtg";

  get<Tags::detail::InputDataSet<Tags::detail::Shift<T>>>(dataset_names_) =
      "/Shift";
  get<Tags::detail::InputDataSet<Tags::detail::Dr<Tags::detail::Shift<T>>>>(
      dataset_names_) = "/DrShift";
  get<Tags::detail::InputDataSet<::Tags::dt<Tags::detail::Shift<T>>>>(
      dataset_names_) = "/DtShift";

  get<Tags::detail::InputDataSet<Tags::detail::Lapse<T>>>(dataset_names_) =
      "/Lapse";
  get<Tags::detail::InputDataSet<Tags::detail::Dr<Tags::detail::Lapse<T>>>>(
      dataset_names_) = "/DrLapse";
  get<Tags::detail::InputDataSet<::Tags::dt<Tags::detail::Lapse<T>>>>(
      dataset_names_) = "/DtLapse";
}

double MetricWorldtubeH5BufferUpdater::update_buffers_for_time(
    const gsl::not_null<Variables<cce_metric_input_tags<ComplexModalVector>>*>
        buffers,
    const gsl::not_null<size_t*> time_span_start,
    const gsl::not_null<size_t*> time_span_end, const double time,
    const size_t computation_l_max, const size_t interpolator_length,
    const size_t buffer_depth) const {
  if (const std::optional<double> next_time =
          common_buffer_checks(time_span_start, time_span_end, time,
                               interpolator_length, buffer_depth);
      next_time.has_value()) {
    return next_time.value();
  }
  // load the desired time spans into the buffers
  // spatial metric
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      tmpl::for_each<Tags::detail::apply_derivs_t<
          Tags::detail::SpatialMetric<ComplexModalVector>>>(
          [this, &i, &j, &buffers, &time_span_start, &time_span_end,
           &computation_l_max](auto tag_v) {
            using tag = typename decltype(tag_v)::type;
            this->update_buffer(
                make_not_null(&get<tag>(*buffers).get(i, j)),
                cce_data_file_.get<h5::Dat>(detail::dataset_name_for_component(
                    get<Tags::detail::InputDataSet<tag>>(dataset_names_), i,
                    j)),
                computation_l_max, *time_span_start, *time_span_end);
            cce_data_file_.close_current_object();
          });
    }
    // shift
    tmpl::for_each<
        Tags::detail::apply_derivs_t<Tags::detail::Shift<ComplexModalVector>>>(
        [this, &i, &buffers, &time_span_start, &time_span_end,
         &computation_l_max](auto tag_v) {
          using tag = typename decltype(tag_v)::type;
          this->update_buffer(
              make_not_null(&get<tag>(*buffers).get(i)),
              cce_data_file_.get<h5::Dat>(detail::dataset_name_for_component(
                  get<Tags::detail::InputDataSet<tag>>(dataset_names_), i)),
              computation_l_max, *time_span_start, *time_span_end);
          cce_data_file_.close_current_object();
        });
  }
  // lapse
  tmpl::for_each<
      Tags::detail::apply_derivs_t<Tags::detail::Lapse<ComplexModalVector>>>(
      [this, &buffers, &time_span_start, &time_span_end,
       &computation_l_max](auto tag_v) {
        using tag = typename decltype(tag_v)::type;
        this->update_buffer(
            make_not_null(&get(get<tag>(*buffers))),
            cce_data_file_.get<h5::Dat>(detail::dataset_name_for_component(
                get<Tags::detail::InputDataSet<tag>>(dataset_names_))),
            computation_l_max, *time_span_start, *time_span_end);
        cce_data_file_.close_current_object();
      });
  // the next time an update will be required
  return time_buffer_[std::min(*time_span_end - interpolator_length + 1,
                               time_buffer_.size() - 1)];
}

double MetricWorldtubeH5BufferUpdater::update_buffers_for_time(
    const gsl::not_null<Variables<cce_metric_input_tags<ComplexDataVector>>*>
        buffers,
    const gsl::not_null<size_t*> time_span_start,
    const gsl::not_null<size_t*> time_span_end, const double time,
    const size_t computation_l_max, const size_t interpolator_length,
    const size_t buffer_depth) const {
  if (const std::optional<double> next_time =
          common_buffer_checks(time_span_start, time_span_end, time,
                               interpolator_length, buffer_depth);
      next_time.has_value()) {
    return next_time.value();
  }
  // load the desired time spans into the buffers
  // spatial metric
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = i; j < 3; ++j) {
      tmpl::for_each<Tags::detail::apply_derivs_t<
          Tags::detail::SpatialMetric<ComplexDataVector>>>(
          [this, &i, &j, &buffers, &time_span_start, &time_span_end,
           &computation_l_max](auto tag_v) {
            using tag = typename decltype(tag_v)::type;
            this->update_buffer(
                make_not_null(&get<tag>(*buffers).get(i, j)),
                cce_data_file_.get<h5::Dat>(detail::dataset_name_for_component(
                    get<Tags::detail::InputDataSet<tag>>(dataset_names_), i,
                    j)),
                computation_l_max, *time_span_start, *time_span_end);
            cce_data_file_.close_current_object();
          });
    }
    // shift
    tmpl::for_each<
        Tags::detail::apply_derivs_t<Tags::detail::Shift<ComplexDataVector>>>(
        [this, &i, &buffers, &time_span_start, &time_span_end,
         &computation_l_max](auto tag_v) {
          using tag = typename decltype(tag_v)::type;
          this->update_buffer(
              make_not_null(&get<tag>(*buffers).get(i)),
              cce_data_file_.get<h5::Dat>(detail::dataset_name_for_component(
                  get<Tags::detail::InputDataSet<tag>>(dataset_names_), i)),
              computation_l_max, *time_span_start, *time_span_end);
          cce_data_file_.close_current_object();
        });
  }
  // lapse
  tmpl::for_each<
      Tags::detail::apply_derivs_t<Tags::detail::Lapse<ComplexDataVector>>>(
      [this, &buffers, &time_span_start, &time_span_end,
       &computation_l_max](auto tag_v) {
        using tag = typename decltype(tag_v)::type;
        this->update_buffer(
            make_not_null(&get(get<tag>(*buffers))),
            cce_data_file_.get<h5::Dat>(detail::dataset_name_for_component(
                get<Tags::detail::InputDataSet<tag>>(dataset_names_))),
            computation_l_max, *time_span_start, *time_span_end);
        cce_data_file_.close_current_object();
      });
  // the next time an update will be required
  return time_buffer_[std::min(*time_span_end - interpolator_length + 1,
                               time_buffer_.size() - 1)];
}

std::unique_ptr<
    WorldtubeBufferUpdater<cce_metric_input_tags<ComplexModalVector>>>
MetricWorldtubeH5BufferUpdater::get_clone() const {
  return std::make_unique<MetricWorldtubeH5BufferUpdater>(
      MetricWorldtubeH5BufferUpdater{filename_});
}

bool MetricWorldtubeH5BufferUpdater::time_is_outside_range(
    const double time) const {
  return time < time_buffer_[0] or time > time_buffer_[time_buffer_.size() - 1];
}

void MetricWorldtubeH5BufferUpdater::pup(PUP::er& p) {
  p | time_buffer_;
  p | has_version_history_;
  p | filename_;
  p | file_is_from_spec_;
  p | l_max_;
  p | extraction_radius_;
  p | is_modal_data_;
  p | dataset_names_;
  if (p.isUnpacking()) {
    cce_data_file_ = h5::H5File<h5::AccessType::ReadOnly>{filename_};
  }
}

std::optional<double> MetricWorldtubeH5BufferUpdater::common_buffer_checks(
    const gsl::not_null<size_t*> time_span_start,
    const gsl::not_null<size_t*> time_span_end, const double time,
    const size_t interpolator_length, const size_t buffer_depth) const {
  if (*time_span_end >= time_buffer_.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (*time_span_end > interpolator_length and
      time_buffer_[*time_span_end - interpolator_length] > time) {
    // the next time an update will be required
    return time_buffer_[*time_span_end - interpolator_length + 1];
  }
  // find the time spans that are needed
  auto new_span_pair = detail::create_span_for_time_value(
      time, buffer_depth, interpolator_length, 0, time_buffer_.size(),
      time_buffer_);
  *time_span_start = new_span_pair.first;
  *time_span_end = new_span_pair.second;

  return std::nullopt;
}

template <typename T>
void MetricWorldtubeH5BufferUpdater::update_buffer(
    const gsl::not_null<T*> buffer_to_update, const h5::Dat& read_data,
    const size_t computation_l_max, const size_t time_span_start,
    const size_t time_span_end) const {
  const size_t number_of_columns = read_data.get_dimensions()[1];
  if (UNLIKELY(buffer_to_update->size() != (time_span_end - time_span_start) *
                                               square(computation_l_max + 1))) {
    ERROR("Incorrect storage size for the data to be loaded in.");
  }
  auto cols = alg::iota(std::vector<size_t>(number_of_columns - 1), 1_st);
  auto data_matrix =
      read_data.get_data_subset<std::vector<std::vector<double>>>(
          cols, time_span_start, time_span_end - time_span_start);

  *buffer_to_update = 0.0;
  for (size_t time_row = 0; time_row < time_span_end - time_span_start;
       ++time_row) {
    // If we have nodal data, we can just copy it into the buffer
    if constexpr (std::is_same_v<T, ComplexDataVector>) {
      ASSERT(not is_modal_data_,
             "Buffer supplied to MetricWorldtubeH5Buffer is for nodal data, "
             "but the class was constructed for modal data.");
      // NOLINTNEXTLINE
      DataVector raw_nodal_data{data_matrix[time_row].data() + 1, cols.size()};
      ComplexDataVector nodal_data_view{
          buffer_to_update->data() + time_row * cols.size(), cols.size()};
      nodal_data_view = std::complex<double>(1.0, 0.0) * raw_nodal_data;

      continue;
    }

    for (int l = 0; l <= static_cast<int>(std::min(computation_l_max, l_max_));
         ++l) {
      for (int m = -l; m <= l; ++m) {
        // -m because SpEC format is stored in decending m.
        const int em = file_is_from_spec_ ? -m : m;
        (*buffer_to_update)[Spectral::Swsh::goldberg_mode_index(
                                computation_l_max, static_cast<size_t>(l), m) *
                                (time_span_end - time_span_start) +
                            time_row] =
            std::complex<double>(
                data_matrix[time_row]
                           [2 * Spectral::Swsh::goldberg_mode_index(
                                    l_max_, static_cast<size_t>(l), em)],
                data_matrix[time_row]
                           [2 * Spectral::Swsh::goldberg_mode_index(
                                    l_max_, static_cast<size_t>(l), em) +
                            1]);
      }
    }
  }
}

BondiWorldtubeH5BufferUpdater::BondiWorldtubeH5BufferUpdater(
    const std::string& cce_data_filename,
    const std::optional<double> extraction_radius, const bool is_modal_data)
    : is_modal_data_(is_modal_data),
      cce_data_file_{cce_data_filename},
      filename_{cce_data_filename} {
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::BondiBeta>>>(dataset_names_) =
      "Beta";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::BondiU>>>(dataset_names_) = "U";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::BondiQ>>>(dataset_names_) = "Q";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::BondiW>>>(dataset_names_) = "W";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::BondiJ>>>(dataset_names_) = "J";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::Dr<Tags::BondiJ>>>>(
      dataset_names_) = "DrJ";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::Du<Tags::BondiJ>>>>(
      dataset_names_) = "H";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::BondiR>>>(dataset_names_) = "R";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::Du<Tags::BondiR>>>>(
      dataset_names_) = "DuR";

  // the extraction radius is typically not used in the Bondi system, so we
  // don't error if it isn't parsed from the filename. Instead, we'll just error
  // if the invalid extraction radius value is ever retrieved using
  // `get_extraction_radius`.
  extraction_radius_ =
      Cce::get_extraction_radius(cce_data_filename, extraction_radius, false);

  detail::set_time_buffer_and_lmax(make_not_null(&time_buffer_), l_max_,
                                   cce_data_file_.get<h5::Dat>("/U"), false,
                                   is_modal_data_);
  cce_data_file_.close_current_object();
}

double BondiWorldtubeH5BufferUpdater::update_buffers_for_time(
    const gsl::not_null<Variables<cce_bondi_input_tags>*> buffers,
    const gsl::not_null<size_t*> time_span_start,
    const gsl::not_null<size_t*> time_span_end, const double time,
    const size_t computation_l_max, const size_t interpolator_length,
    const size_t buffer_depth) const {
  return detail::update_buffers_for_time<cce_bondi_input_tags>(
      buffers, time_span_start, time_span_end, time, computation_l_max, l_max_,
      interpolator_length, buffer_depth, time_buffer_, dataset_names_,
      cce_data_file_, is_modal_data_);
}

void BondiWorldtubeH5BufferUpdater::pup(PUP::er& p) {
  p | time_buffer_;
  p | filename_;
  p | l_max_;
  p | extraction_radius_;
  p | is_modal_data_;
  p | dataset_names_;
  if (p.isUnpacking()) {
    cce_data_file_ = h5::H5File<h5::AccessType::ReadOnly>{filename_};
  }
}

KleinGordonWorldtubeH5BufferUpdater::KleinGordonWorldtubeH5BufferUpdater(
    const std::string& cce_data_filename,
    const std::optional<double> extraction_radius)
    : cce_data_file_{cce_data_filename}, filename_{cce_data_filename} {
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::KleinGordonPsi>>>(
      dataset_names_) = "KGPsi";
  get<Tags::detail::InputDataSet<
      Spectral::Swsh::Tags::SwshTransform<Tags::KleinGordonPi>>>(
      dataset_names_) = "dtKGPsi";

  // the extraction radius is typically not used in the Klein-Gordon system, so
  // we don't error if it isn't parsed from the filename. Instead, we'll just
  // error if the invalid extraction radius value is ever retrieved using
  // `get_extraction_radius`.
  extraction_radius_ =
      Cce::get_extraction_radius(cce_data_filename, extraction_radius, false);

  detail::set_time_buffer_and_lmax(make_not_null(&time_buffer_), l_max_,
                                   cce_data_file_.get<h5::Dat>("/KGPsi"), true,
                                   true);
  cce_data_file_.close_current_object();
}

double KleinGordonWorldtubeH5BufferUpdater::update_buffers_for_time(
    const gsl::not_null<Variables<klein_gordon_input_tags>*> buffers,
    const gsl::not_null<size_t*> time_span_start,
    const gsl::not_null<size_t*> time_span_end, const double time,
    const size_t computation_l_max, const size_t interpolator_length,
    const size_t buffer_depth) const {
  return detail::update_buffers_for_time<klein_gordon_input_tags>(
      buffers, time_span_start, time_span_end, time, computation_l_max, l_max_,
      interpolator_length, buffer_depth, time_buffer_, dataset_names_,
      cce_data_file_, true);
}

void KleinGordonWorldtubeH5BufferUpdater::pup(PUP::er& p) {
  p | time_buffer_;
  p | filename_;
  p | l_max_;
  p | extraction_radius_;
  p | dataset_names_;
  if (p.isUnpacking()) {
    cce_data_file_ = h5::H5File<h5::AccessType::ReadOnly>{filename_};
  }
}

PUP::able::PUP_ID MetricWorldtubeH5BufferUpdater::my_PUP_ID = 0;
PUP::able::PUP_ID BondiWorldtubeH5BufferUpdater::my_PUP_ID = 0;
PUP::able::PUP_ID KleinGordonWorldtubeH5BufferUpdater::my_PUP_ID = 0;  // NOLINT

#define SPIN(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                       \
  template void detail::update_buffer_with_modal_data<SPIN(data)>( \
      const gsl::not_null<ComplexModalVector*> buffer_to_update,   \
      const h5::Dat& read_data, const size_t computation_l_max,    \
      const size_t l_max, const size_t time_span_start,            \
      const size_t time_span_end, const bool is_modal_data);

GENERATE_INSTANTIATIONS(INSTANTIATE, (0, 1, 2))

#undef INSTANTIATE
#undef SPIN

#define VEC_TYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                 \
  template void MetricWorldtubeH5BufferUpdater::set_names<VEC_TYPE(data)>(); \
  template void MetricWorldtubeH5BufferUpdater::update_buffer(               \
      gsl::not_null<VEC_TYPE(data)*> buffer_to_update,                       \
      const h5::Dat& read_data, size_t computation_l_max,                    \
      size_t time_span_start, size_t time_span_end) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (ComplexModalVector, ComplexDataVector))

#undef INSTANTIATE
#undef VEC_TYPE
}  // namespace Cce
