// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <unordered_map>
#include <vector>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Storage.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
namespace domain {
class BlockId;
}  // namespace domain
template <size_t VolumeDim>
class ElementId;
template <typename IdType, typename DataType>
class IdPair;
/// \endcond

namespace ah {
// Assumes the vars in `interpolated_vars_storage` is already properly sized.
bool interpolate_to_target_points(
    gsl::not_null<Storage::InterpolatedVars*> interpolated_vars_storage,
    const std::unordered_map<ElementId<3>, Storage::VolumeVars>&
        element_ids_and_volume_vars,
    const std::vector<std::optional<
        IdPair<domain::BlockId, tnsr::I<double, 3, ::Frame::BlockLogical>>>>&
        block_logical_coords);
}  // namespace ah
