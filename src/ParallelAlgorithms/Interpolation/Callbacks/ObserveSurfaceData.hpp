// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/YlmSpherepack.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/PostInterpolationCallback.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp {
namespace callbacks {

/// \brief post_interpolation_callback that outputs
/// 2D "volume" data on a surface.
///
/// Uses:
/// - Metavariables
///   - `temporal_id`
/// - DataBox:
///   - `TagsToObserve` (each tag must be a Scalar<DataVector>)
///
/// Conforms to the intrp::protocols::PostInterpolationCallback protocol
///
/// For requirements on InterpolationTargetTag, see
/// intrp::protocols::InterpolationTargetTag
template <typename TagsToObserve, typename InterpolationTargetTag,
          typename HorizonFrame>
struct ObserveSurfaceData
    : tt::ConformsTo<intrp::protocols::PostInterpolationCallback> {
  static constexpr double fill_invalid_points_with =
      std::numeric_limits<double>::quiet_NaN();

  using const_global_cache_tags = tmpl::list<observers::Tags::SurfaceFileName>;

  template <typename DbTags, typename Metavariables, typename TemporalId>
  static void apply(const db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const TemporalId& temporal_id) {
    const Strahlkorper<HorizonFrame>& strahlkorper =
        get<StrahlkorperTags::Strahlkorper<HorizonFrame>>(box);
    const YlmSpherepack& ylm = strahlkorper.ylm_spherepack();
    const auto& inertial_strahlkorper_coords =
        get<StrahlkorperTags::CartesianCoords<::Frame::Inertial>>(box);
    const auto& grid_strahlkorper_coords =
        get<StrahlkorperTags::CartesianCoords<::Frame::Grid>>(box);

    // Output the inertial-frame coordinates of the Stralhlkorper.
    // Note that these coordinates are not
    // Spherepack-evenly-distributed over the inertial-frame sphere
    // (they are Spherepack-evenly-distributed over the HorizonFrame
    // sphere).
    std::vector<TensorComponent> tensor_components{
        {"InertialCoordinates_x"s, get<0>(inertial_strahlkorper_coords)},
        {"InertialCoordinates_y"s, get<1>(inertial_strahlkorper_coords)},
        {"InertialCoordinates_z"s, get<2>(inertial_strahlkorper_coords)},
        {"GridCoordinates_x"s, get<0>(grid_strahlkorper_coords)},
        {"GridCoordinates_y"s, get<1>(grid_strahlkorper_coords)},
        {"GridCoordinates_z"s, get<2>(grid_strahlkorper_coords)}};

    // Output each tag if it is a scalar. Otherwise, throw a compile-time
    // error. This could be generalized to handle tensors of nonzero rank by
    // looping over the components, so each component could be visualized
    // separately as a scalar. But in practice, this generalization is
    // probably unnecessary, because Strahlkorpers are typically only
    // visualized with scalar quantities (used set the color at different
    // points on the surface).
    tmpl::for_each<TagsToObserve>([&box, &tensor_components](auto tag_v) {
      using Tag = tmpl::type_from<decltype(tag_v)>;
      static_assert(std::is_same_v<typename Tag::type, Scalar<DataVector>>,
                    "Each tag in TagsToObserve must be a Scalar<DataVector>. "
                    "This could be generalized if needed; see the code comment "
                    "above the static_assert.");
      tensor_components.push_back({db::tag_name<Tag>(), get(get<Tag>(box))});
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
    const observers::ObservationId& observation_id = observers::ObservationId(
        InterpolationTarget_detail::get_temporal_id_value(temporal_id),
        subfile_path + ".vol");

    auto& proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);

    // We call this on proxy[0] because the 0th element of a NodeGroup is
    // always guaranteed to be present.
    Parallel::threaded_action<observers::ThreadedActions::WriteVolumeData>(
        proxy[0], Parallel::get<observers::Tags::SurfaceFileName>(cache),
        subfile_path, observation_id,
        std::vector<ElementVolumeData>{{extents_vector, tensor_components,
                                        bases_vector, quadratures_vector,
                                        surface_name}});
  }
};
}  // namespace callbacks
}  // namespace intrp
