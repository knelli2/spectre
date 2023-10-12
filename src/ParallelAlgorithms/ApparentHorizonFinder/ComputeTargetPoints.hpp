// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "DataStructures/IdPair.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

/// \cond
class DataVector;
class FastFlow;
namespace domain::Tags {
template <size_t Dim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace ylm {
template <typename Frame>
class Strahlkorper;
}  // namespace ylm
namespace Frame {
struct NoFrame;
struct BlockLogical;
struct Grid;
struct Distorted;
struct Inertial;
}  // namespace Frame
/// \endcond

namespace ah {
tnsr::I<DataVector, 3, Frame::NoFrame> compute_target_points(
    const ylm::Strahlkorper<Frame::NoFrame>& old_horizon,
    const ::FastFlow& fast_flow);

namespace detail {
template <typename Frame>
tnsr::I<DataVector, 3, Frame> points_in_frame(
    const tnsr::I<DataVector, 3, Frame::NoFrame>& points_no_frame) {
  tnsr::I<DataVector, 3, Frame> result{};
  for (size_t i = 0; i < 3; i++) {
    result.get(i) = DataVector{points_no_frame.get(i).data(),
                               points_no_frame.get(i).size()};
  }

  return result;
}
}  // namespace detail

// Use the output from `compute_target_points`
template <typename Metavariables>
std::vector<std::optional<
    IdPair<domain::BlockId, tnsr::I<double, 3, ::Frame::BlockLogical>>>>
compute_block_logical_target_points(
    const Parallel::GlobalCache<Metavariables>& cache, const double time,
    const tnsr::I<DataVector, 3, Frame::NoFrame>& points_no_frame,
    const std::string& frame) {
  const Domain<3>& domain = get<domain::Tags::Domain<3>>(cache);

  if (frame == "Grid") {
    // Frame is grid frame, so don't need any FunctionsOfTime,
    // whether or not the maps are time_dependent.
    return ::block_logical_coordinates(
        domain, detail::points_in_frame<Frame::Grid>(points_no_frame));
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
        return ::block_logical_coordinates(
            domain, detail::points_in_frame<Frame::Distorted>(points_no_frame),
            time, functions_of_time);
      } else {
        return ::block_logical_coordinates(
            domain, detail::points_in_frame<Frame::Inertial>(points_no_frame),
            time, functions_of_time);
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
  return ::block_logical_coordinates(
      domain, detail::points_in_frame<Frame::Grid>(points_no_frame));
}
}  // namespace ah
