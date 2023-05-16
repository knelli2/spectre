// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/Cce/IsBoundaryElement.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits.hpp"

#include "Parallel/Printf.hpp"

/// \cond
template <class Metavariables, class PhaseDepActionList>
struct DgElementArray;
namespace logging::Tags {
template <typename OptionsGroup>
struct Verbosity;
}
/// \endcond

namespace Cce {
namespace Actions {

template <typename CceComponent, bool FirstTime>
struct SendNextTimeToCcm {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    if (not Cce::is_outer_boundary_element<false>(make_not_null(&inboxes), box,
                                                  cache)) {
      // if (Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(
      //         cache) >= ::Verbosity::Debug) {
      //   Parallel::printf(
      //       "SendNextTimeToCcm %s: Not a boundary element; not sending "
      //       "anything.\n",
      //       array_index);
      // }
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    const TimeStepId& next_time =
        db::get<::Tags::Next<::Tags::TimeStepId>>(box);
    const TimeStepId& current_time = db::get<::Tags::TimeStepId>(box);

    if (Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(cache) >=
        ::Verbosity::Debug) {
      Parallel::printf(
          "SendNextTimeToCcm %s: Sending next time %f at current time %f\n",
          array_index, next_time.substep_time(), current_time.substep_time());
    }

    Parallel::receive_data<Cce::ReceiveTags::GhNextTimeToCcm>(
        Parallel::get_parallel_component<CceComponent>(cache), current_time,
        // QUESTION: true??????
        std::make_pair(array_index, next_time), false);

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Actions
}  // namespace Cce
