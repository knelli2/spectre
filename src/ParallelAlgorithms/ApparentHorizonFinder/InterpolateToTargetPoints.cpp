// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/ApparentHorizonFinder/Actions/TryToInterpolate.hpp"

#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/IdPair.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/ElementLogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Storage.hpp"
#include "Utilities/Gsl.hpp"

namespace ah {
void interpolate_to_target_points(
    const gsl::not_null<Storage::InterpolatedVars*> interpolated_vars_storage,
    const std::unordered_map<ElementId<3>, Storage::VolumeVars>&
        element_ids_and_volume_vars,
    const std::vector<std::optional<
        IdPair<domain::BlockId, tnsr::I<double, 3, ::Frame::BlockLogical>>>>&
        block_logical_coords) {
  // Get list of ElementIds that have not yet been interpolated.
  std::vector<ElementId<3>> element_ids;
  element_ids.reserve(element_ids_and_volume_vars.size());

  for (const auto& [element_id, volume_vars] : element_ids_and_volume_vars) {
    (void)volume_vars;
    // Have we interpolated this element before?
    if (interpolated_vars_storage->interpolation_is_done_for_these_elements
            .find(element_id) ==
        interpolated_vars_storage->interpolation_is_done_for_these_elements
            .end()) {
      interpolated_vars_storage->interpolation_is_done_for_these_elements
          .emplace(element_id);
      element_ids.emplace_back(element_id);
    }
  }

  // Get element logical coordinates.
  const std::unordered_map<ElementId<3>, ElementLogicalCoordHolder<3>>
      element_coord_holders =
          element_logical_coordinates(element_ids, block_logical_coords);

  auto& indices_interpolated_to_thus_far =
      interpolated_vars->indicies_interpolated_to_thus_far;
  auto& interpolated_vars = interpolated_vars_storage->vars;
  const size_t number_of_grid_points =
      interpolated_vars.number_of_grid_points();
  const size_t number_of_components =
      interpolated_vars.number_of_independent_components;

  // Loop over every element that can be interpolated to and do the
  // interpolation
  for (const auto& [element_id, element_coord_holder] : element_coord_holders) {
    // Create the interpolant
    intrp::Irregular<3> interpolator(
        volume_info.mesh, element_coord_holder.element_logical_coords);

    const Storage::VolumeVars& volume_vars =
        element_ids_and_volume_vars.at(element_id);

    const auto& global_offsets = element_coord_holder.offsets;
    auto interpolated_vars_on_single_element =
        interpolator.interpolate(volume_vars.volume_vars_for_horizon_finding);

    // Loop over points in this element and fill their values in the big buffer
    const size_t number_of_interpolated_points = global_offsets.size();
    for (size_t i = 0; i < number_of_interpolated_points; ++i) {
      // If a point is on the boundary of two (or more) elements, it is possible
      // that we have received data for this point from a different element.
      // This will rarely occur, but it does occur, e.g. when a point is exactly
      // on some symmetry boundary (such as the x-y plane) and this symmetry
      // boundary is exactly the boundary between two elements.  If this
      // happens, we accept the first duplicated point, and we ignore subsequent
      // duplicated points.  The points are easy to keep track
      // of because global_offsets uniquely identifies them.
      if (indices_interpolated_to_thus_far.insert(global_offsets[i]).second) {
        for (size_t c = 0; c < number_of_components; ++c) {
          // clang-tidy: no pointer arithmetic
          // NOLINTBEGIN
          interpolated_vars
              .data()[global_offsets[i] + c * number_of_grid_points] =
              interpolated_vars_on_single_element
                  .data()[i + c * number_of_interpolated_points];
          // NOLINTEND
        }
      }
    }

    // Add this element to the finished ones
    interpolated_vars_storage->interpolation_is_done_for_these_elements.insert(
        element_id);
  }

  return interpolated_vars->indicies_interpolated_to_thus_far.size() ==
         block_logical_coords.size();
}
}  // namespace ah
