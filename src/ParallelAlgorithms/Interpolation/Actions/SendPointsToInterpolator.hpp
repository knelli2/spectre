// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/Logging/Tags.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp {
namespace Actions {
/// \ingroup ActionsGroup
/// \brief Sets up points on an `InterpolationTarget` at a new `temporal_id`
/// and sends these points to an `Interpolator`.
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
  template <
      typename ParallelComponent, typename DbTags, typename Metavariables,
      typename ArrayIndex, typename TemporalId,
      Requires<tmpl::list_contains_v<DbTags, Tags::TemporalIds<TemporalId>>> =
          nullptr>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const TemporalId& temporal_id) {
    const auto& verbosity =
        db::get<logging::Tags::Verbosity<InterpolationTargetTag>>(box);
    auto coords = InterpolationTarget_detail::block_logical_coords<
        InterpolationTargetTag>(box, cache, temporal_id);
    InterpolationTarget_detail::set_up_interpolation<InterpolationTargetTag>(
        make_not_null(&box), temporal_id, coords);
    auto& receiver_proxy =
        Parallel::get_parallel_component<Interpolator<Metavariables>>(cache);
    if (verbosity > ::Verbosity::Verbose) {
      Parallel::printf(
          "%s: t=%.6g: SendPointsToInterpolator: "
          "Calling Actions::ReceivePoints and exiting\n",
          pretty_type::short_name<InterpolationTargetTag>(),
          InterpolationTarget_detail::get_temporal_id_value(temporal_id));
    }
    Parallel::simple_action<Actions::ReceivePoints<InterpolationTargetTag>>(
        receiver_proxy, temporal_id, std::move(coords));
  }
};

}  // namespace Actions
}  // namespace intrp
