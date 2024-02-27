// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_set>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/Algorithms/AlgortihmArray.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Actions/InitializeInterpolationTarget.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits.hpp"

/// \cond
namespace intrp {
namespace Actions {
template <typename Metavariables, typename InterpolationTargetTag>
struct InitializeInterpolationTarget;
}  // namespace Actions
}  // namespace intrp
namespace evolution {
struct EventsAndDenseTriggers;
namespace Tags {
struct EventsAndDenseTriggers;
}  // namespace Tags
namespace OptionTags {
struct EventsAndDenseTriggers;
}  // namespace OptionTags
}  // namespace evolution
/// \endcond

namespace intrp2 {
template <class Metavariables>
struct InterpolationTarget {
  using metavariables = Metavariables;
  static constexpr bool include_dense_triggers =
      Metavariables::intrp::include_dense_triggers;
  using chare_type = ::Parallel::Algorithms::Array;
  using array_index = std::string;

  using const_global_cache_tags = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<intrp2::Actions::InitializeInterpolationTarget>>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  using array_allocation_tags = tmpl::conditional_t<
      include_dense_triggers,
      tmpl::list<evolution::OptionTags::EventsAndDenseTriggers>, tmpl::list<>>;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const std::unordered_set<size_t>& procs_to_ignore = {}) {
    allocate_array_impl(global_cache, initialization_items,
                        evolution::EventsAndDenseTriggers{}, procs_to_ignore);
  }

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const evolution::EventsAndDenseTriggers& events_and_dense_triggers,
      const std::unordered_set<size_t>& procs_to_ignore = {}) {
    allocate_array_impl(global_cache, initialization_items,
                        events_and_dense_triggers, procs_to_ignore);
  }

  static void execute_next_phase(
      Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<InterpolationTarget2>(local_cache)
        .start_phase(next_phase);
  };

 private:
  static void allocate_array_impl(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const evolution::EventsAndDenseTriggers& events_and_dense_triggers,
      const std::unordered_set<size_t>& procs_to_ignore);
};

template <class Metavariables>
void InterpolationTarget2<Metavariables>::allocate_array_impl(
    Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
    const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
        initialization_items,
    const evolution::EventsAndDenseTriggers& events_and_dense_triggers,
    const std::unordered_set<size_t>& procs_to_ignore) {
  auto& local_cache = *Parallel::local_branch(global_cache);
  auto& interpolation_target =
      Parallel::get_parallel_component<InterpolationTarget2>(local_cache);

  const EventsAndTriggers& events_and_triggers =
      Parallel::get<::Tags::EventsAndTriggers>(local_cache);

  std::unordered_set<std::string> interpolation_events{};

  const auto is_interpolation_event =
      [&interpolation_events](const ::Event& event) {
        if (event.is_interpolation()) {
          interpolation_events.insert(event.name());
        }
      };

  events_and_triggers.for_each_event(is_interpolation_event);

  if constexpr (include_dense_triggers) {
    events_and_dense_triggers.for_each_event(is_interpolation_event);
  } else {
    (void)events_and_dense_triggers;
  }

  const size_t number_of_procs = Parallel::number_of_procs<size_t>(local_cache);

  // TODO: This needs a better method, but this probably requires rewriting the
  // singleton allocation part of the ResourceInfo
  size_t which_proc = 0;
  for (const std::string& array_component_name : interpolation_events) {
    while (procs_to_ignore.find(which_proc) != procs_to_ignore.end()) {
      which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
    }
    interpolation_target(array_component_name)
        .insert(global_cache, initialization_items, which_proc);

    which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
  }
}
}  // namespace intrp2
