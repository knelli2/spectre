// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <vector>

template <size_t Dim>
class ElementId;

namespace domain {
/// \brief Computes the Z-curve index of a given `ElementId`
///
/// \details The Z-curve index is computed by interleaving the bits of the
/// `ElementId`'s `Segment` indices. Here is a sketch of a 2D block with 4x2
/// elements, with bit indices and the resulting z-curve:
///
/// \code
///        x-->
///        00  01  10  11
/// y  0 |  0   2   4   6
/// |    |
/// v  1 |  1   3   5   7
/// \endcode
///
/// \param element_id the `ElementId` for which to compute the Z-curve index
template <size_t Dim>
size_t z_curve_index_from_element_id(const ElementId<Dim>& element_id);

/// \brief Computes the `Segment` indices of an `Element` from a given Z-curve
/// index and block refinements
///
/// \details The `Segment` indices are computed by extracting the bits of the
/// `Segment` indices from the Z-curve index. Here is a sketch of a 2D block
/// with 4x2 elements, with bit indices and the resulting z-curve:
///
/// \code
///        x-->
///        00  01  10  11
/// y  0 |  0   2   4   6
/// |    |
/// v  1 |  1   3   5   7
/// \endcode
///
/// \param z_curve_index the Z-curve index for which to compute the `Segment`
/// indices
/// \param initial_ref_levs the refinements of the block that the `ElementId`
/// belongs to
template <size_t Dim>
std::array<size_t, Dim> segment_indices_from_z_curve_index(
    const size_t z_curve_index,
    const std::array<size_t, Dim>& initial_ref_levs);

/// \brief Create the `ElementId`s of a single `Block` ordered by their Z-curve
/// index
///
/// \details For details on `ElementId`'s Z-curve indices, see
/// `domain::segment_indices_from_z_curve_index`
///
/// \param block_id the `Block` number
/// \param initial_ref_levs the refinements levels of the `Block`
/// \param grid_index the grid index of the `Block`
template <size_t Dim>
std::vector<ElementId<Dim>> initial_element_ids_in_z_curve_order(
    size_t block_id, std::array<size_t, Dim> initial_ref_levs,
    size_t grid_index = 0);
}  // namespace domain
