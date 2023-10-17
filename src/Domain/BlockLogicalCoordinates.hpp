// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/IdPair.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Parallel/GlobalCache.hpp"

/// \cond
template <size_t VolumeDim>
class Domain;
template <size_t VolumeDim>
class Block;
namespace domain::Tags {
template <size_t VolumeDim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

template <size_t Dim>
using BlockLogicalCoords = std::optional<
    IdPair<domain::BlockId, tnsr::I<double, Dim, Frame::BlockLogical>>>;

/// @{
/// \ingroup ComputationalDomainGroup
///
/// Computes the block logical coordinates and the containing `BlockId` of
/// a set of points, given coordinates in a particular frame.
///
/// \details Returns a std::vector<std::optional<IdPair<BlockId,coords>>>,
/// where the vector runs over the points and is indexed in the same order as
/// the input coordinates `x`. For each point, the `IdPair` holds the
/// block logical coords of that point and the `BlockId` of the `Block` that
/// contains that point.
/// The std::optional is invalid if the point is not in any Block.
/// If a point is on a shared boundary of two or more `Block`s, it is
/// returned only once, and is considered to belong to the `Block`
/// with the smaller `BlockId`.
///
/// The `block_logical_coordinates_single_point` function will search the passed
/// in block for the passed in coordinate and return the logical coordinates of
/// that point. It will return a `std::nullopt` if it can't find the point in
/// that block.
///
/// \warning Since map inverses can involve numerical roundoff error, care must
/// be taken with points on shared block boundaries. They will be assigned to
/// the first block (by block ID) that contains the point _within roundoff
/// error_. Therefore, be advised to use the logical coordinates returned by
/// this function, which are guaranteed to be in [-1, 1] and can be safely
/// passed along to `element_logical_coordinates`.
///
/// \warning `block_logical_coordinates` with x in
/// `::Frame::Distorted` ignores all `Block`s that lack a distorted
/// frame, and it will return std::nullopt for points that lie outside
/// all distorted-frame-endowed `Block`s. This is what is expected for
/// typical use cases.  This means that `block_logical_coordinates`
/// does not assume that grid and distorted frames are equal in
/// `Block`s that lack a distorted frame.
template <size_t Dim, typename Fr>
auto block_logical_coordinates(
    const Domain<Dim>& domain, const tnsr::I<DataVector, Dim, Fr>& x,
    double time = std::numeric_limits<double>::signaling_NaN(),
    const domain::FunctionsOfTimeMap& functions_of_time = {})
    -> std::vector<BlockLogicalCoords<Dim>>;

template <size_t Dim, typename Fr>
std::optional<tnsr::I<double, Dim, ::Frame::BlockLogical>>
block_logical_coordinates_single_point(
    const tnsr::I<double, Dim, Fr>& input_point, const Block<Dim>& block,
    double time = std::numeric_limits<double>::signaling_NaN(),
    const domain::FunctionsOfTimeMap& functions_of_time = {});
/// @}

namespace block_logical_detail {
template <size_t Dim, typename Frame>
void points_in_frame(
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame>*> points_with_frame,
    const tnsr::I<DataVector, Dim, ::Frame::NoFrame>& points_no_frame) {
  for (size_t i = 0; i < Dim; i++) {
    points_with_frame->get(i).set_data_ref(
        const_cast<DataVector&>(points_no_frame.get(i)).data(),
        points_no_frame.get(i).size());
  }
}
}  // namespace block_logical_detail

template <size_t Dim, typename Metavariables>
std::vector<std::optional<
    IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>
block_logical_coordinates_in_frame(
    const Parallel::GlobalCache<Metavariables>& cache, const double time,
    const tnsr::I<DataVector, Dim, Frame::NoFrame>& points_no_frame,
    const std::string& frame) {
  const Domain<Dim>& domain = get<domain::Tags::Domain<Dim>>(cache);

  if (frame == "Grid") {
    // Frame is grid frame, so don't need any FunctionsOfTime,
    // whether or not the maps are time_dependent.
    tnsr::I<DataVector, Dim, Frame::Grid> points_with_frame{};
    block_logical_detail::points_in_frame(make_not_null(&points_with_frame),
                                          points_no_frame);
    return ::block_logical_coordinates(domain, points_with_frame);
  }

  if (domain.is_time_dependent()) {
    if constexpr (Parallel::is_in_mutable_global_cache<
                      Metavariables, domain::Tags::FunctionsOfTime>) {
      // Whoever calls this function when the maps are time-dependent is
      // responsible for ensuring that functions_of_time are up to date at the
      // time.
      const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
      // It's either Distorted or Inertial now
      if (frame == "Distorted") {
        tnsr::I<DataVector, Dim, Frame::Distorted> points_with_frame{};
        block_logical_detail::points_in_frame(make_not_null(&points_with_frame),
                                              points_no_frame);
        return ::block_logical_coordinates(domain, points_with_frame, time,
                                           functions_of_time);
      } else {
        tnsr::I<DataVector, Dim, Frame::Inertial> points_with_frame{};
        block_logical_detail::points_in_frame(make_not_null(&points_with_frame),
                                              points_no_frame);
        return ::block_logical_coordinates(domain, points_with_frame, time,
                                           functions_of_time);
      }
    } else {
      // We error here because the maps are time-dependent, yet
      // the cache does not contain FunctionsOfTime.  It would be
      // nice to make this a compile-time error; however, we want
      // the code to compile for the completely time-independent
      // case where there are no FunctionsOfTime in the cache at
      // all.  Unfortunately, checking whether the maps are
      // time-dependent is currently not constexpr.
      ERROR(
          "There is a time-dependent CoordinateMap in at least one "
          "of the Blocks, but FunctionsOfTime are not in the "
          "GlobalCache.  If you intend to use a time-dependent "
          "CoordinateMap, please add FunctionsOfTime to the GlobalCache.");
    }
  }

  // Time-independent case. Regardless of what the frame actually is, Grid,
  // Distorted, and Inertial are all the same so we just choose Grid.
  tnsr::I<DataVector, Dim, Frame::Grid> points_with_frame{};
  block_logical_detail::points_in_frame(make_not_null(&points_with_frame),
                                        points_no_frame);
  return ::block_logical_coordinates(domain, points_with_frame);
}
