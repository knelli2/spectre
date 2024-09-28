// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/RadialTranslation.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <ostream>
#include <pup.h>
#include <pup_stl.h>
#include <type_traits>
#include <unordered_set>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/DereferenceWrapper.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain::CoordinateMaps::TimeDependent {

template <size_t Dim>
RadialTranslation<Dim>::RadialTranslation(std::string function_of_time_name,
                                          const std::array<double, Dim>& center)
    : inner_f_of_t_name_(std::move(function_of_time_name)),
      outer_f_of_t_name_(std::move(function_of_time_name)),
      functions_of_time_equal_(true),
      f_of_t_names_({inner_f_of_t_name_}),
      center_(center) {}

template <size_t Dim>
RadialTranslation<Dim>::RadialTranslation(
    std::string inner_function_of_time_name,
    std::string outer_function_of_time_name, double inner_radius,
    double outer_radius, const std::array<double, Dim>& center)
    : inner_f_of_t_name_(std::move(inner_function_of_time_name)),
      outer_f_of_t_name_(std::move(outer_function_of_time_name)),
      functions_of_time_equal_(inner_f_of_t_name_ == outer_f_of_t_name_),
      f_of_t_names_({inner_f_of_t_name_, outer_f_of_t_name_}),
      inner_radius_(inner_radius),
      outer_radius_(outer_radius),
      center_(center) {
  if (outer_radius_ <= inner_radius_) {
    ERROR("Outer radius " << outer_radius_
                          << " must be larger than inner radius "
                          << inner_radius_);
  }
  one_over_radii_difference_ = 1.0 / (outer_radius_ - inner_radius_);
}

template <size_t Dim>
template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, Dim> RadialTranslation<Dim>::operator()(
    const std::array<T, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  return function_and_velocity_helper(source_coords, time, functions_of_time,
                                      0);
}

template <size_t Dim>
std::optional<std::array<double, Dim>> RadialTranslation<Dim>::inverse(
    const std::array<double, Dim>& target_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  std::array<double, Dim> centered_target_coords = target_coords - center_;
  const double centered_target_radius = magnitude(centered_target_coords);
  double scaling_factor = 0.0;

  if (functions_of_time_equal_) {
    const DataVector inner_function_of_time =
        functions_of_time.at(inner_f_of_t_name_)->func(time)[0];
    ASSERT(inner_function_of_time.size() == 1,
           "The dimension of the function of time ("
               << inner_function_of_time.size() << ") must be 1.");

    scaling_factor = centered_target_radius - inner_function_of_time[0];

  } else {
    const DataVector inner_function_of_time =
        functions_of_time.at(inner_f_of_t_name_)->func(time)[0];
    ASSERT(inner_function_of_time.size() == 1,
           "The dimension of the inner function of time ("
               << inner_function_of_time.size() << ") must be 1.");
    const DataVector outer_function_of_time =
        functions_of_time.at(outer_f_of_t_name_)->func(time)[0];
    ASSERT(outer_function_of_time.size() == 1,
           "The dimension of the outer function of time ("
               << outer_function_of_time.size() << ") must be 1.");

    scaling_factor =
        (centered_target_radius -
         one_over_radii_difference_ *
             (outer_radius_ * inner_function_of_time[0] -
              inner_radius_ * outer_function_of_time[0])) /
        (1.0 + one_over_radii_difference_ *
                   (outer_function_of_time[0] - inner_function_of_time[0]));
  }

  const double eps = std::numeric_limits<double>::epsilon() * 100.0;
  const double source_radius = scaling_factor * centered_target_radius;
  if (source_radius + eps < inner_radius_ or
      source_radius - eps > outer_radius_) {
    return std::nullopt;
  }

  return scaling_factor * centered_target_coords / centered_target_radius +
         center_;
}

template <size_t Dim>
template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, Dim>
RadialTranslation<Dim>::frame_velocity(
    const std::array<T, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  return function_and_velocity_helper(source_coords, time, functions_of_time,
                                      1);
}

template <size_t Dim>
template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>
RadialTranslation<Dim>::jacobian(
    const std::array<T, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  using CVT = tt::remove_cvref_wrap_t<T>;
  auto result = make_with_value<tnsr::Ij<CVT, Dim, Frame::NoFrame>>(
      dereference_wrapper(source_coords[0]), 0.0);
  std::array<CVT, Dim> centered_coords{};
  for (size_t i = 0; i < Dim; i++) {
    gsl::at(centered_coords, i) =
        gsl::at(source_coords, i) - gsl::at(center_, i);
  }
  CVT one_over_centered_radius = magnitude(centered_coords);
  for (size_t i = 0; i < get_size(one_over_centered_radius); i++) {
    ASSERT(get_element(one_over_centered_radius, i) > 0.0,
           "Source coord is at the center " << center_
                                            << "! This should not happen.");
  }
  one_over_centered_radius = 1.0 / one_over_centered_radius;
  const CVT one_over_centered_radius_cubed = cube(one_over_centered_radius);

  const DataVector inner_function_of_time =
      functions_of_time.at(inner_f_of_t_name_)->func(time)[0];
  ASSERT(inner_function_of_time.size() == 1,
         "The dimension of the inner function of time ("
             << inner_function_of_time.size() << ") must be 1.");

  CVT factor_for_cross_terms{};
  CVT factor_for_delta{};

  if (functions_of_time_equal_) {
    factor_for_cross_terms =
        inner_function_of_time[0] * one_over_centered_radius_cubed;
    factor_for_delta =
        1.0 + inner_function_of_time[0] * one_over_centered_radius;
  } else {
    const DataVector outer_function_of_time =
        functions_of_time.at(outer_f_of_t_name_)->func(time)[0];
    ASSERT(outer_function_of_time.size() == 1,
           "The dimension of the outer function of time ("
               << outer_function_of_time.size() << ") must be 1.");
    check_bounding_radii(inner_function_of_time, outer_function_of_time);

    factor_for_cross_terms = one_over_radii_difference_ *
                             (outer_radius_ * inner_function_of_time[0] -
                              inner_radius_ * outer_function_of_time[0]) *
                             one_over_centered_radius_cubed;
    factor_for_delta =
        1.0 +
        one_over_radii_difference_ *
            (outer_function_of_time[0] - inner_function_of_time[0]) +
        factor_for_cross_terms * one_over_centered_radius;
  }

  for (size_t i = 0; i < Dim; i++) {
    result.get(i, i) += factor_for_delta;
    for (size_t j = 0; j < Dim; j++) {
      result.get(i, j) -= factor_for_cross_terms * gsl::at(centered_coords, i) *
                          gsl::at(centered_coords, j);
    }
  }

  return result;
}
template <size_t Dim>
template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, Dim, Frame::NoFrame>
RadialTranslation<Dim>::inv_jacobian(
    const std::array<T, Dim>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const {
  return determinant_and_inverse(
             jacobian(source_coords, time, functions_of_time))
      .second;
}

template <size_t Dim>
template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, Dim>
RadialTranslation<Dim>::function_and_velocity_helper(
    const std::array<T, Dim>& source_coords, double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    const size_t function_or_deriv_index) const {
  using CVT = tt::remove_cvref_wrap_t<T>;
  const auto inner_func_and_deriv_of_time =
      functions_of_time.at(inner_f_of_t_name_)->func_and_deriv(time);
  ASSERT(inner_func_and_deriv_of_time[0].size() == 1,
         "The dimension of the inner function of time ("
             << inner_func_and_deriv_of_time[0].size() << ") must be 1.");

  // TODO: Have points inside inner radius and outside outer radius just have
  // rigid radial translations.

  std::array<CVT, Dim> centered_source_coords{};
  for (size_t i = 0; i < Dim; i++) {
    gsl::at(centered_source_coords, i) =
        gsl::at(source_coords, i) - gsl::at(center_, i);
  }

  const CVT centered_radius = magnitude(centered_source_coords);
  CVT scaling_factor{};

  if (functions_of_time_equal_) {
    scaling_factor = centered_radius + gsl::at(inner_func_and_deriv_of_time,
                                               function_or_deriv_index)[0];
  } else {
    const auto outer_func_and_deriv_of_time =
        functions_of_time.at(outer_f_of_t_name_)->func_and_deriv(time);
    ASSERT(outer_func_and_deriv_of_time[0].size() == 1,
           "The dimension of the outer function of time ("
               << outer_func_and_deriv_of_time[0].size() << ") must be 1.");
    check_bounding_radii(inner_func_and_deriv_of_time[0],
                         outer_func_and_deriv_of_time[0]);

    scaling_factor =
        one_over_radii_difference_ *
        (outer_radius_ *
             gsl::at(inner_func_and_deriv_of_time, function_or_deriv_index)[0] -
         inner_radius_ *
             gsl::at(outer_func_and_deriv_of_time, function_or_deriv_index)[0] +
         (gsl::at(outer_func_and_deriv_of_time, function_or_deriv_index)[0] -
          gsl::at(inner_func_and_deriv_of_time, function_or_deriv_index)[0]) *
             centered_radius);
    if (function_or_deriv_index == 0) {
      scaling_factor += centered_radius;
    }
  }

  return scaling_factor * centered_source_coords / centered_radius + center_;
}

template <size_t Dim>
void RadialTranslation<Dim>::check_bounding_radii(
    const DataVector& inner_function_of_time,
    const DataVector& outer_function_of_time) const {
  if (inner_radius_ + inner_function_of_time[0] + 1.0e-8 >
      outer_radius_ + outer_function_of_time[0]) {
    ERROR(
        "Inner radius has been mapped beyond the outer radius. This is a bug "
        "and should not happen.");
  }
}

template <size_t Dim>
void RadialTranslation<Dim>::pup(PUP::er& p) {
  size_t version = 0;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  p | inner_f_of_t_name_;
  p | outer_f_of_t_name_;
  p | functions_of_time_equal_;
  p | f_of_t_names_;
  p | inner_radius_;
  p | outer_radius_;
  p | one_over_radii_difference_;
  p | center_;
}

template <size_t Dim>
bool operator==(const RadialTranslation<Dim>& lhs,
                const RadialTranslation<Dim>& rhs) {
  return lhs.inner_f_of_t_name_ == rhs.inner_f_of_t_name_ and
         lhs.outer_f_of_t_name_ == rhs.outer_f_of_t_name_ and
         lhs.center_ == rhs.center_ and
         lhs.functions_of_time_equal_ == rhs.functions_of_time_equal_ and
         // Avoids comparing NaNs
         (not lhs.functions_of_time_equal_ and
          lhs.inner_radius_ == rhs.inner_radius_ and
          lhs.outer_radius_ == rhs.outer_radius_);
}

// Explicit instantiations
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                    \
  template class RadialTranslation<DIM(data)>;                  \
  template bool operator==(const RadialTranslation<DIM(data)>&, \
                           const RadialTranslation<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef INSTANTIATE

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                                                \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)>      \
  RadialTranslation<DIM(data)>::operator()(                                 \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const;                                         \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)>      \
  RadialTranslation<DIM(data)>::frame_velocity(                             \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const;                                         \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),        \
                    Frame::NoFrame>                                         \
  RadialTranslation<DIM(data)>::jacobian(                                   \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const;                                         \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),        \
                    Frame::NoFrame>                                         \
  RadialTranslation<DIM(data)>::inv_jacobian(                               \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (2, 3),
                        (double, DataVector,
                         std::reference_wrapper<const double>,
                         std::reference_wrapper<const DataVector>))
#undef DIM
#undef DTYPE
#undef INSTANTIATE
}  // namespace domain::CoordinateMaps::TimeDependent
