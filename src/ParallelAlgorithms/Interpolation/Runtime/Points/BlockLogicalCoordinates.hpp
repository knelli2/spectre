// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/IdPair.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TypeTraits/CreateHasTypeAlias.hpp"

namespace intrp2 {
template <size_t Dim>
using BlockCoords = std::vector<std::optional<
    IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>;

namespace detail {
CREATE_HAS_TYPE_ALIAS(argument_tags)
CREATE_HAS_TYPE_ALIAS_V(argument_tags)
}  // namespace detail

/*!
 * \brief Given a struct `Points` that conforms to the
 * `intrp2::protocols::Points` protocol, get the block logical coordinates at a
 * specific time.
 *
 * \param box Either a DataBox or ObservationBox
 * \param cache The GlobalCache
 * \param time Time of this interpolation
 * \param frame Frame that this interpolation will be in
 * \param points A class that conforms to the `intrp2::protocols::Points`
 * protocol
 * \return std::optional<BlockCoords> The output of the
 * `::block_logical_coordinates()` function.
 */
template <typename BoxType, typename Metavariables, typename Points>
std::optional<BlockCoords<Metavariables::volume_dim>> block_logical_coordinates(
    const BoxType& box, const Parallel::GlobalCache<Metavariables>& cache,
    const double time, const std::string& frame, const Points& points) {
  if constexpr (detail::has_argument_tags_v<Points>) {
    return db::apply(points, box, cache, time, frame);
  } else {
    (void)box;
    return ::block_logical_coordinates_in_frame(
        cache, time, points.target_points_no_frame(), frame);
  }
}
}  // namespace intrp2
