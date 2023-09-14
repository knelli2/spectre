// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <functional>
#include <optional>
#include <unordered_set>
#include <vector>

#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/ResourceInfo.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/OptionTags.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace ah {
/*!
 * \brief Horizon Finder parallel component
 */
template <class Metavariables>
struct HorizonFinderComponent {
  static constexpr size_t volume_dim = Metavariables::volume_dim;

  using chare_type = Parallel::Algorithms::Array;
  using metavariables = Metavariables;
  using array_index = int;

  using const_global_cache_tags = tmpl::list<>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization, tmpl::list<>,
                             Parallel::Actions::TerminatePhase>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const std::unordered_set<size_t>& procs_to_ignore = {});

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<HorizonFinderComponent>(local_cache)
        .start_phase(next_phase);
  }
};

template <class Metavariables>
void HorizonFinderComponent<Metavariables>::allocate_array(
    Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
    const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
        initialization_items,
    const std::unordered_set<size_t>& /*procs_to_ignore*/) {
  auto& local_cache = *Parallel::local_branch(global_cache);
  auto& horizon_finder =
      Parallel::get_parallel_component<HorizonFinderComponent>(local_cache);

  //   const Parallel::ResourceInfo<Metavariables>& resource_info =
  //       local_cache.get_resource_info();

  const std::optional<std::vector<HorizonFinderOptions>>&
      horizons_from_options =
          tuples::get<Tags::HorizonFinders>(initialization_items);

  if (horizons_from_options.has_value()) {
    for (const HorizonFinderOptions& options : horizons_from_options.value()) {
      // FIXME: Don't always insert on core 0. Use ResourceInfo? But this isn't
      // a singleton anymore
      horizon_finder(static_cast<int>(options.object_label))
          .insert(global_cache, initialization_items, 0);
    }
  }

  horizon_finder.doneInserting();
}
}  // namespace ah
