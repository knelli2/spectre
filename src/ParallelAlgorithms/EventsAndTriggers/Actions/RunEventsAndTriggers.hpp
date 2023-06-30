// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace Tags {
struct TimeStepId;
}  // namespace Tags
/// \endcond

namespace Actions {
/// \ingroup ActionsGroup
/// \ingroup EventsAndTriggersGroup
/// \brief Run the events and triggers
///
/// Uses:
/// - GlobalCache: the EventsAndTriggers tag, as required by
///   events and triggers
/// - DataBox: as required by events and triggers
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies: nothing
struct RunEventsAndTriggers {
  using const_global_cache_tags = tmpl::list<Tags::EventsAndTriggers>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const component) {
    // Checking the slab status in this way will work for
    // non-evolutions and does not require linking against any
    // evolution-related libraries because everything is duck-typed
    // off the DataBox.
    if constexpr (db::tag_is_retrievable_v<Tags::TimeStepId,
                                           db::DataBox<DbTags>>) {
      if (not db::get<Tags::TimeStepId>(box).is_at_slab_boundary()) {
        return {Parallel::AlgorithmExecution::Continue, std::nullopt};
      }
    }

    Parallel::get<Tags::EventsAndTriggers>(cache).run_events(
        box, cache, array_index, component);

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Actions