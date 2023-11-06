// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

#include "DataStructures/IdPair.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Utilities/StdHelpers.hpp"

namespace intrp {
namespace Actions {
/// \ingroup ActionsGroup
/// \brief Sets up points on an `InterpolationTarget` at a new `temporal_id`
/// and sends these points to an `Interpolator`.
///
/// The `iteration` parameter tags each set of points so the `Interpolator`
/// knows which are newer points and which are older points.
///
/// \see `intrp::Actions::ReceivePoints`
///
/// Uses:
/// - DataBox:
///   - `domain::Tags::Domain<3>`
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies:
///   - `Tags::IndicesOfFilledInterpPoints`
///   - `Tags::IndicesOfInvalidInterpPoints`
///   - `Tags::InterpolatedVars<InterpolationTargetTag, TemporalId>`
///
/// For requirements on InterpolationTargetTag, see InterpolationTarget
template <typename InterpolationTargetTag>
struct SendPointsToInterpolator {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex, typename TemporalId>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const TemporalId& temporal_id,
                    const size_t iteration = 0_st) {
    auto coords = InterpolationTarget_detail::block_logical_coords<
        InterpolationTargetTag>(box, cache, temporal_id);
    for (size_t i = 0; i < coords.size(); i++) {
      const auto& block_coord = coords[i];
      if (not block_coord.has_value()) {
        using ::operator<<;
        const std::string preface =
            InterpolationTarget_detail::target_output_prefix<
                SendPointsToInterpolator, InterpolationTargetTag>(temporal_id);
        const auto distorted_points =
            InterpolationTargetTag::compute_target_points::points(
                box, tmpl::type_<Metavariables>{}, temporal_id);
        const std::array distorted_point{
            get<0>(distorted_points)[i],
            get<1>(distorted_points)[i],
            get<2>(distorted_points)[i],
        };
        ERROR_NO_TRACE(preface << ", Not all block logical coordinates exist.\n"
                                  " Index of first coord that doesn't exist = "
                               << i
                               << "\n Distorted point = " << distorted_point
                               << "\n Full vector: " << coords
                               << "\n Distorted points: " << distorted_points);
      }
    }
    InterpolationTarget_detail::set_up_interpolation<InterpolationTargetTag>(
        make_not_null(&box), temporal_id, coords);
    auto& receiver_proxy =
        Parallel::get_parallel_component<Interpolator<Metavariables>>(cache);
    Parallel::simple_action<Actions::ReceivePoints<InterpolationTargetTag>>(
        receiver_proxy, temporal_id, std::move(coords), iteration);
    if (Parallel::get<intrp::Tags::Verbosity>(cache) >= ::Verbosity::Debug) {
      Parallel::printf(
          "%s, Sending points to interpolator.\n",
          InterpolationTarget_detail::target_output_prefix<
              SendPointsToInterpolator, InterpolationTargetTag>(temporal_id));
    }
  }
};

}  // namespace Actions
}  // namespace intrp
