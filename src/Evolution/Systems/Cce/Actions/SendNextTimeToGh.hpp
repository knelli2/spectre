// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

#include "Parallel/Printf.hpp"

/// \cond
template <class Metavariables, class PhaseDepActionList>
struct DgElementArray;
struct DebugToggle;
/// \endcond

namespace Cce {
namespace Actions {

template <typename CceComponent, bool FirstTime>
struct SendNextTimeToGh {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const TimeStepId& next_time =
        db::get<::Tags::Next<::Tags::TimeStepId>>(box);
    const TimeStepId& current_time = db::get<::Tags::TimeStepId>(box);

    using gh_dg_element_array =
        tmpl::filter<typename Metavariables::component_list,
                     tt::is_a<::DgElementArray, tmpl::_1>>;

    if constexpr (tmpl::size<gh_dg_element_array>::value == 1) {
      if (Parallel::get<DebugToggle>(cache)) {
        Parallel::printf(
            "SendNextTimeToGh: Sending next time %f at current time %f\n",
            next_time.substep_time().value(),
            current_time.substep_time().value());
      }

      // if constexpr (FirstTime) {
      //   Parallel::receive_data<Cce::ReceiveTags::CcmNextTimeToGH>(
      //     Parallel::get_parallel_component<tmpl::front<gh_dg_element_array>>(
      //           cache),
      //       TimeStepId{current_time.time_runs_forward(),
      //                  current_time.slab_number(),
      //                  Time{current_time.substep_time().slab().retreat(),
      //                       current_time.substep_time().fraction()}},
      //       current_time, false);
      // }

      if constexpr (FirstTime) {
        Parallel::receive_data<Cce::ReceiveTags::CcmNextTimeToGH>(
            Parallel::get_parallel_component<tmpl::front<gh_dg_element_array>>(
                cache),
            current_time, current_time, false);
      } else {
        Parallel::receive_data<Cce::ReceiveTags::CcmNextTimeToGH>(
            Parallel::get_parallel_component<tmpl::front<gh_dg_element_array>>(
                cache),
            current_time, next_time, false);
      }
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Actions
}  // namespace Cce
