// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/Algorithms/AlgorithmSingletonDeclarations.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/MemoryMonitor/Initialize.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup ParallelGroup
 * \brief Singleton parallel component used for monitoring memory usage of other
 * parallel components
 */
template <class Metavariables>
struct MemoryMonitor {
  using chare_type = Parallel::Algorithms::Singleton;

  using metavariables = Metavariables;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename metavariables::Phase, metavariables::Phase::Initialization,
      tmpl::list<
          Actions::SetupDataBox,
          Initialization::Actions::InitializeMemoryMonitor<Metavariables>,
          Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;

  using initialization_tags = Parallel::get_initialization_tags<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const typename Metavariables::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<MemoryMonitor<Metavariables>>(local_cache)
        .start_phase(next_phase);
  }
};
