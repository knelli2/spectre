// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataBox/ValidateSelection.hpp"
#include "Domain/StrahlkorperTransformations.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
struct Inertial;
struct Grid;
struct Distorted;
}  // namespace Frame
namespace intrp2::callbacks {
template <typename TagsToObserve>
struct ObserveDataOnStrahlkorper;
}  // namespace intrp2::callbacks
/// \endcond

namespace intrp2 {
namespace callbacks {
namespace detail {
// Fills the legend and row of spherical harmonic data to write to disk
//
// The number of coefficients to write is based on `max_l`, the maximum value
// that the input `strahlkorper` could possibly have. When
// `strahlkorper.l_max() < max_l`, coefficients with \f$l\f$ higher than
// `strahlkorper.l_max()` will simply be zero. Assuming the same `max_l` is
// always used for a given surface, we will always write the same number of
// columns for each row, as `max_l` sets the number of columns to write
template <typename Frame>
void fill_ylm_legend_and_data(gsl::not_null<std::vector<std::string>*> legend,
                              gsl::not_null<std::vector<double>*> data,
                              const ylm::Strahlkorper<Frame>& strahlkorper,
                              double time, size_t max_l);
}  // namespace detail

/// \brief post_interpolation_callback that outputs 2D "volume" data on a
/// surface and the surface's spherical harmonic data
///
/// \details
/// Uses:
/// - Metavariables
///   - `temporal_id`
/// - DataBox:
///   - `TagsToObserve` (each tag must be a Scalar<DataVector>)
///
/// Conforms to the intrp2::protocols::PostInterpolationCallback protocol
///
/// For requirements on InterpolationTargetTag, see
/// intrp2::protocols::InterpolationTargetTag
///
/// The columns of spherical harmonic data written take the form
///
/// \code
/// [Time, {Frame}ExpansionCenter_x, {Frame}ExpansionCenter_y,
/// {Frame}ExpansionCenter_z, Lmax, coef(0,0), ... coef(Lmax,Lmax)]
/// \endcode
///
/// where `coef(l,m)` refers to the strahlkorper coefficients stored and defined
/// by `ylm::Strahlkorper::coefficients() const`. It is assumed that
/// \f$l_{max} = m_{max}\f$.
///
/// \note Currently, \f$l_{max}\f$ for a given surface does not change over the
/// course of the simulation, which means that the total number of columns of
/// coefficients that we need to write is also constant. The current
/// implementation of writing the coefficients at one time assumes \f$l_{max}\f$
/// of a surface remains constant. If and when in the future functionality for
/// an adaptive \f$l_{max}\f$ is added, the implementation for writing the
/// coefficients will need to be updated to account for this. One possible way
/// to address this is to have a known maximum \f$l_{max}\f$ for a given surface
/// and write all coefficients up to that maximum \f$l_{max}\f$.
template <typename Target, typename... TagsToObserve>
struct ObserveDataOnStrahlkorper<tmpl::list<TagsToObserve...>>
    : public Callback<Target>, protocols::Callback {
  using tags_to_observe_on_target = tmpl::list<>;
  using non_observation_tags_on_target = tmpl::list<>;
  using volume_compute_tags = tmpl::list<>;

  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help = {
        "The name of the subfile inside the HDF5 file without an extension and "
        "without a preceding '/'."};
  };
  struct ValuesToObserve {
    using type = std::vector<std::string>;
    static constexpr Options::String help = {
        "List specifying the name of each value to observe."};
  };

  using options = tmpl::list<SubfileName, ValuesToObserve>;

  using available_tags_to_observe = tmpl::list<TagsToObserve...>;

  static constexpr Options::String help = {
      "Observe values (doubles) over a whole surface. These are usually "
      "reduced or integrated quantities."};

  ObserveDataOnStrahlkorper(std::string subfile_name,
                            const std::vector<std::string>& values_to_observe,
                            const Options::Context& context = {})
      : subfile_name_(std::move(subfile_name)) {
    db::validate_selection<available_tags_to_observe>(values_to_observe);
    for (const auto& value : values_to_observe) {
      values_to_observe_.insert(value);
    }
  }

  void pup(PUP::er& p) override {
    p | subfile_name_;
    p | values_to_observe_;
  }

  template <typename Metavariables>
  static void apply(const db::Access& access,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const double time) {
    const ylm::Strahlkorper<Frame::NoFrame>& strahlkorper =
        get<ylm::Tags::Strahlkorper<Frame::NoFrame>>(access);
    const ylm::Spherepack& ylm = strahlkorper.ylm_spherepack();

    std::vector<TensorComponent> tensor_components;

    const std::string& frame = db::get<intrp2::Tags::Frame>(access);
    const auto& coords_in_frame = db::get<intrp2::Tags::Points<3>>(access);
    tensor_components.push_back(
        {frame + "Coordinates_x"s, get<0>(coords_in_frame)});
    tensor_components.push_back(
        {frame + "Coordinates_y"s, get<1>(coords_in_frame)});
    tensor_components.push_back(
        {frame + "Coordinates_z"s, get<2>(coords_in_frame)});

    if (frame != "Inertial") {
      tnsr::I<DataVector, 3, Frame::Inertial> inertial_coords{};

      const auto& domain = Parallel::get<domain::Tags::Domain<3>>(cache);
      const auto& functions_of_time =
          Parallel::get<domain::Tags::FunctionOfTime>(cache);

      if (frame == "Distorted") {
        ylm::Strahlkorper<Frame::Distorted> strahlkorper_dist_frame{
            strahlkorper.l_max(), strahlkorper.radius(),
            strahlkorper.expansion_center()};
        strahlkorper_dist_frame.coefficients().set_data_ref(
            strahlkorper.coefficients().data(),
            strahlkorper.coefficients().size());
        strahlkorper_coords_in_different_frame(make_not_null(&inertial_coords),
                                               strahlkorper_dist_frame, domain,
                                               functions_of_time, time);
      } else {
        ASSERT(frame == "Grid", "Unknown frame '"
                                    << frame
                                    << "' for ObserveDataOnStrahlkorper.");
        ylm::Strahlkorper<Frame::Grid> strahlkorper_grid_frame{
            strahlkorper.l_max(), strahlkorper.radius(),
            strahlkorper.expansion_center()};
        strahlkorper_grid_frame.coefficients().set_data_ref(
            strahlkorper.coefficients().data(),
            strahlkorper.coefficients().size());
        strahlkorper_coords_in_different_frame(make_not_null(&inertial_coords),
                                               strahlkorper_grid_frame, domain,
                                               functions_of_time, time);
      }

      tensor_components.push_back(
          {"InertialCoordinates_x"s, get<0>(inertial_coords)});
      tensor_components.push_back(
          {"InertialCoordinates_y"s, get<1>(inertial_coords)});
      tensor_components.push_back(
          {"InertialCoordinates_z"s, get<2>(inertial_coords)});
    }

    const std::unordered_set<std::string>& vars_to_observe =
        db::get<intrp2::Tags::VarsToObserve>(access);

    // Output each tag if it is a scalar. Otherwise, throw a compile-time
    // error. This could be generalized to handle tensors of nonzero rank by
    // looping over the components, so each component could be visualized
    // separately as a scalar. But in practice, this generalization is
    // probably unnecessary, because Strahlkorpers are typically only
    // visualized with scalar quantities (used set the color at different
    // points on the surface).
    tmpl::for_each<TagsToObserve>(
        [&box, &tensor_components, &vars_to_observe](auto tag_v) {
          using Tag = tmpl::type_from<decltype(tag_v)>;
          const auto tag_name = db::tag_name<Tag>();
          if (vars_to_observe.count(tag_name) != 1) {
            return;
          }

          const auto& tensor = get<Tag>(box);
          for (size_t i = 0; i < tensor.size(); ++i) {
            tensor_components.emplace_back(
                tag_name + tensor.component_suffix(i), tensor[i]);
          }
        });

    const std::string& surface_name =
        pretty_type::name<InterpolationTargetTag>();
    const std::string subfile_path{std::string{"/"} + surface_name};
    const std::vector<size_t> extents_vector{
        {ylm.physical_extents()[0], ylm.physical_extents()[1]}};
    const std::vector<Spectral::Basis> bases_vector{
        2, Spectral::Basis::SphericalHarmonic};
    const std::vector<Spectral::Quadrature> quadratures_vector{
        {Spectral::Quadrature::Gauss, Spectral::Quadrature::Equiangular}};
    const observers::ObservationId observation_id{time, subfile_path + ".vol"};

    auto& proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);

    // We call this on proxy[0] because the 0th element of a NodeGroup is
    // always guaranteed to be present.
    Parallel::threaded_action<observers::ThreadedActions::WriteVolumeData>(
        proxy[0], Parallel::get<observers::Tags::SurfaceFileName>(cache),
        subfile_path, observation_id,
        std::vector<ElementVolumeData>{{surface_name, tensor_components,
                                        extents_vector, bases_vector,
                                        quadratures_vector}});

    std::vector<std::string> ylm_legend{};
    std::vector<double> ylm_data{};
    // The number of coefficients written will be (l_max + 1)^2 where l_max is
    // the current value of l_max for this surface's strahlkorper. Because l_max
    // remains constant, the number of coefficient columns written does, too. In
    // the future when l_max is adaptive, instead of passing in the current
    // l_max of the strahlkorper, we could pass in the maximum value that l_max
    // could be to ensure that we (a) have enough columns to write all the
    // coefficients regardless of the current value of l_max and (b) write a
    // constant number of columns for each row of data regardless of the current
    // l_max.
    detail::fill_ylm_legend_and_data(make_not_null(&ylm_legend),
                                     make_not_null(&ylm_data), strahlkorper,
                                     time, strahlkorper.l_max());

    const std::string ylm_subfile_name{std::string{"/"} + surface_name +
                                       "_Ylm"};

    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        proxy[0], ylm_subfile_name, std::move(ylm_legend),
        std::make_tuple(std::move(ylm_data)));
  }

  const std::unordered_set<std::string>& observables() const override {
    return values_to_observe_;
  }

 private:
  std::string subfile_name_{};
  std::unordered_set<std::string> values_to_observe_{};
};
}  // namespace callbacks
}  // namespace intrp2
