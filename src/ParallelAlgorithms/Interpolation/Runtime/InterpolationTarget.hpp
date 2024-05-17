// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_set>

#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Protocols/ArrayElementsAllocator.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Metavariables.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace intrp2 {
struct EventsAndTriggersAllocator
    : public tt::ConformsTo<Parallel::protocols::ArrayElementsAllocator> {
  template <typename ParallelComponent>
  using array_allocation_tags = tmpl::list<>;

  template <typename ParallelComponent, typename Metavariables,
            typename... InitializationTags>
  static void apply(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::TaggedTuple<InitializationTags...>& initialization_items,
      const tuples::tagged_tuple_from_typelist<
          typename ParallelComponent::array_allocation_tags>&
      /*array_allocation_items*/
      = {},
      const std::unordered_set<size_t>& procs_to_ignore = {}) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    auto& interpolation_target =
        Parallel::get_parallel_component<ParallelComponent>(local_cache);

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

    const size_t number_of_procs =
        Parallel::number_of_procs<size_t>(local_cache);

    // TODO: This needs a better method, but this probably requires rewriting
    // the singleton allocation part of the ResourceInfo
    size_t which_proc = 0;
    for (const std::string& array_component_name : interpolation_events) {
      while (procs_to_ignore.contains(which_proc)) {
        which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
      }
      interpolation_target(array_component_name)
          .insert(global_cache, initialization_items, which_proc);

      which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
    }
  }
};

template <class Metavariables>
struct InterpolationTargets {
  using metavariables = Metavariables;
  using IntrpMetavars = typename Metavariables::intrp;
  static_assert(tt::assert_conforms_to_v<IntrpMetavars,
                                         intrp2::protocols::Metavariables>);
  using ElementsAllocator = typename IntrpMetavars::elements_allocator;
  static_assert(
      tt::assert_conforms_to_v<ElementsAllocator,
                               Parallel::protocols::ArrayElementsAllocator>);
  using chare_type = ::Parallel::Algorithms::Array;
  using array_index = std::string;

  using const_global_cache_tags = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<intrp2::Actions::InitializeInterpolationTarget>>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  using array_allocation_tags =
      typename ElementsAllocator::template array_allocation_tags<
          InterpolationTargets>;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const tuples::tagged_tuple_from_typelist<array_allocation_tags>&
          array_allocation_items = {},
      const std::unordered_set<size_t>& procs_to_ignore = {}) {
    ElementsAllocator::template apply<InterpolationTargets>(
        global_cache, initialization_items, array_allocation_items,
        procs_to_ignore);
  }

  static void execute_next_phase(
      Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<InterpolationTargets>(local_cache)
        .start_phase(next_phase);
  };
};
}  // namespace intrp2
