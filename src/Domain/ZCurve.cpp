// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/ZCurve.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"

namespace domain {

template <size_t Dim>
size_t z_curve_index_from_element_id(const ElementId<Dim>& element_id) {
  // for the bit manipulation of the element index, we need to interleave the
  // indices in each dimension in order according to how many bits are in the
  // index representation. This variable stores the refinement level and
  // dimension index in ascending order of refinement level, representing a
  // permutation of the dimensions
  // pair<refinement level, dim index> in order of ascending refinement
  std::array<std::pair<size_t, size_t>, Dim>
      dimension_by_highest_refinement_level;
  for (size_t i = 0; i < Dim; ++i) {
    dimension_by_highest_refinement_level.at(i) =
        std::make_pair(element_id.segment_id(i).refinement_level(), i);
  }
  alg::sort(dimension_by_highest_refinement_level,
            [](const std::pair<size_t, size_t>& lhs,
               const std::pair<size_t, size_t>& rhs) {
              return lhs.first < rhs.first;
            });

  size_t element_order_index = 0;

  // 'gap' the lowest refinement direction bits as:
  // ... x1 x0 -> ... x1 0 0 x0,
  // then bitwise or in 'gap'ed and shifted next-lowest refinement direction
  // bits as:
  // ... y2 y1 y0 -> ... y2 0 y1 x1 0 y0 x0
  // then bitwise or in 'gap'ed and shifted highest-refinement direction bits
  // as:
  // ... z3 z2 z1 z0 -> z3 z2 y2 z1 y1 x1 z0 y0 x0
  // note that we must skip refinement-level 0 dimensions as though they are
  // not present
  size_t leading_gap = 0;
  for (size_t i = 0; i < Dim; ++i) {
    const size_t id_to_gap_and_shift =
        element_id
            .segment_id(
                gsl::at(dimension_by_highest_refinement_level, i).second)
            .index();
    size_t total_gap = leading_gap;
    if (gsl::at(dimension_by_highest_refinement_level, i).first > 0) {
      ++leading_gap;
    }
    for (size_t bit_index = 0;
         bit_index < gsl::at(dimension_by_highest_refinement_level, i).first;
         ++bit_index) {
      // This operation will not overflow for our present use of `ElementId`s.
      // This technique densely assigns an ElementID a unique size_t identifier
      // determining the Morton curve order, and `ElementId` supports refinement
      // levels such that a global index within a block will fit in a 64-bit
      // unsigned integer.
      element_order_index |=
          ((id_to_gap_and_shift & two_to_the(bit_index)) << total_gap);
      for (size_t j = 0; j < Dim; ++j) {
        if (i != j and
            bit_index + 1 <
                gsl::at(dimension_by_highest_refinement_level, j).first) {
          ++total_gap;
        }
      }
    }
  }
  return element_order_index;
}

template <size_t Dim>
std::array<size_t, Dim> segment_indices_from_z_curve_index(
    const size_t z_curve_index,
    const std::array<size_t, Dim>& initial_ref_levs) {
  std::array<std::pair<size_t, size_t>, Dim>
      dimension_by_highest_refinement_level;
  for (size_t i = 0; i < Dim; ++i) {
    dimension_by_highest_refinement_level.at(i) =
        std::make_pair(gsl::at(initial_ref_levs, i), i);
  }
  // dimensions in order of ascending refinement
  alg::sort(dimension_by_highest_refinement_level,
            [](const std::pair<size_t, size_t>& lhs,
               const std::pair<size_t, size_t>& rhs) {
              return lhs.first < rhs.first;
            });

  // result segment indices to compute
  std::array<std::pair<size_t, size_t>, Dim>
      segment_indices_by_highest_refinement_level;
  for (size_t i = 0; i < Dim; ++i) {
    segment_indices_by_highest_refinement_level.at(i) = std::make_pair(
        0, gsl::at(dimension_by_highest_refinement_level, i).second);
  }

  // index of lowest-refined dimension to loop over for extracting ElementId
  // bits (any dim below this index has already had all bits extracted)
  size_t starting_dim_index = 0;
  // current bit position of result segment indices, i.e. how many bits of the
  // segment indices that we've already extracted
  size_t bit_index = 0;
  // Z-curve index value for extracting the result ElementId
  size_t element_order_index = z_curve_index;

  // Extract all but the highest bits of the highest refined dimension that are
  // not shared by the dimension with the next-highest refinement. For example,
  // if we have refinement 2 in x, 3 in y, and 5 in z, this will extract all
  // bits of x and y, but only the first 3 bits of z since the next-highest
  // refinement (y) is 3. i.e. extract only the z2 y2 z1 y1 x1 z0 y0 x0 bits
  // of the whole Z-curve index, z4 z3 z2 y2 z1 y1 x1 z0 y0 x0.
  while (starting_dim_index < Dim - 1) {
    // refinement level of lowest dim we're still extracting bits for
    const size_t refinement_level =
        gsl::at(dimension_by_highest_refinement_level, starting_dim_index)
            .first;
    // number of bits to extract for each dimension we're still extracting bits
    // for
    const size_t num_bits_to_extract = refinement_level - bit_index;
    for (size_t i = 0; i < num_bits_to_extract; i++) {
      for (size_t dim_index = starting_dim_index; dim_index < Dim;
           dim_index++) {
        // extract lowest bit of what remains of the Z-curve index value, shift
        // the bit up to the bit position we're on, and add it to the result
        // segment index
        segment_indices_by_highest_refinement_level.at(dim_index).first |=
            ((element_order_index & 1) << bit_index);
        // clear the lowest bit that we just completed extracting
        element_order_index >>= 1;
      }
      bit_index++;
    }

    starting_dim_index++;
    // move the starting dim index to the dim with the next-highest refinement
    // since we've extracted all bits up to this refinement level for all dims
    while (starting_dim_index < Dim and
           refinement_level == gsl::at(dimension_by_highest_refinement_level,
                                       starting_dim_index)
                                   .first) {
      starting_dim_index++;
    }
  }
  // Extract the remaining highest bits of the highest refined dimension that
  // are not shared by the dimension with the next-highest refinement. For
  // example, if we have refinement 2 in x, 3 in y, and 5 in z, then all bits of
  // x and y and only the first 3 bits of z have already been extracted, so this
  // will just extract the remaining 2 bits for z. i.e. extract only the
  // remaining z4 z3 bits of the whole Z-curve index,
  // z4 z3 z2 y2 z1 y1 x1 z0 y0 x0.
  //
  // If the highest refinement is shared by more than one dimension, this will
  // have no effect since all bits have already been extracted in the loop
  // above.
  segment_indices_by_highest_refinement_level.at(Dim - 1).first |=
      (element_order_index << bit_index);

  // result segment indices reordered with the dimension ordering of the input,
  // i.e. undo the initial sorting of dimension by highest refinement
  std::array<size_t, Dim> result_segment_indices{};
  for (size_t i = 0; i < Dim; i++) {
    gsl::at(result_segment_indices,
            gsl::at(segment_indices_by_highest_refinement_level, i).second) =
        gsl::at(segment_indices_by_highest_refinement_level, i).first;
  }

  return result_segment_indices;
}

template <>
std::vector<ElementId<1>> initial_element_ids_in_z_curve_order<1>(
    const size_t block_id, const std::array<size_t, 1> initial_ref_levs,
    const size_t grid_index) {
  std::vector<ElementId<1>> ids;
  const size_t num_elements = two_to_the(initial_ref_levs[0]);
  ids.reserve(num_elements);

  for (size_t i = 0; i < num_elements; i++) {
    std::array<size_t, 1> segment_indices =
        domain::segment_indices_from_z_curve_index(i, initial_ref_levs);
    SegmentId x_segment_id(initial_ref_levs[0], segment_indices[0]);
    ids.emplace_back(block_id, make_array<1>(x_segment_id), grid_index);
  }

  return ids;
}

template <>
std::vector<ElementId<2>> initial_element_ids_in_z_curve_order<2>(
    const size_t block_id, const std::array<size_t, 2> initial_ref_levs,
    const size_t grid_index) {
  std::vector<ElementId<2>> ids;
  const size_t num_elements =
      two_to_the(initial_ref_levs[0]) * two_to_the(initial_ref_levs[1]);
  ids.reserve(num_elements);

  for (size_t i = 0; i < num_elements; i++) {
    std::array<size_t, 2> segment_indices =
        domain::segment_indices_from_z_curve_index(i, initial_ref_levs);
    SegmentId x_segment_id(initial_ref_levs[0], segment_indices[0]);
    SegmentId y_segment_id(initial_ref_levs[1], segment_indices[1]);
    ids.emplace_back(block_id, make_array(x_segment_id, y_segment_id),
                     grid_index);
  }

  return ids;
}

template <>
std::vector<ElementId<3>> initial_element_ids_in_z_curve_order<3>(
    const size_t block_id, const std::array<size_t, 3> initial_ref_levs,
    const size_t grid_index) {
  std::vector<ElementId<3>> ids;
  const size_t num_elements = two_to_the(initial_ref_levs[0]) *
                              two_to_the(initial_ref_levs[1]) *
                              two_to_the(initial_ref_levs[2]);
  ids.reserve(num_elements);

  for (size_t i = 0; i < num_elements; i++) {
    std::array<size_t, 3> segment_indices =
        domain::segment_indices_from_z_curve_index(i, initial_ref_levs);
    SegmentId x_segment_id(initial_ref_levs[0], segment_indices[0]);
    SegmentId y_segment_id(initial_ref_levs[1], segment_indices[1]);
    SegmentId z_segment_id(initial_ref_levs[2], segment_indices[2]);
    ids.emplace_back(block_id,
                     make_array(x_segment_id, y_segment_id, z_segment_id),
                     grid_index);
  }

  return ids;
}

#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                     \
  template size_t z_curve_index_from_element_id(   \
      const ElementId<GET_DIM(data)>& element_id); \
  template std::array<size_t, GET_DIM(data)>       \
  segment_indices_from_z_curve_index(              \
      const size_t z_curve_index,                  \
      const std::array<size_t, GET_DIM(data)>& initial_ref_levs);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

}  // namespace domain
