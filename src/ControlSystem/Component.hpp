// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <tuple>

#include "ControlSystem/Actions/Initialization.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "IO/Observer/Actions/RegisterSingleton.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/TMPL.hpp"

// This is the component for a _single_ control system.
template <class Metavariables, typename ControlSystem>
struct ControlComponent {
  using chare_type = Parallel::Algorithms::Singleton;

  using metavariables = Metavariables;

  //clang-format off
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename metavariables::Phase, metavariables::Phase::Initialization,
      tmpl::list<Actions::SetupDataBox,
                 Initialization::Actions::InitializeControlSystem<
                     Metavariables, ControlSystem>,
                 Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;

  using initialization_tags = Parallel::get_initialization_tags<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const typename Metavariables::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) noexcept {
    auto& local_cache = *(global_cache.ckLocalBranch());
    Parallel::get_parallel_component<
        ControlComponent<Metavariables, ControlSystem>>(local_cache)
        .start_phase(next_phase);
  }
};

namespace control_system {
/// \ingroup ControlSystemGroup
/// List of control componenets to be added to the component list of the
/// metavars
template <typename Metavariables, typename ControlSystems>
using control_components = tmpl::transform<
    ControlSystems,
    tmpl::bind<ControlComponent, tmpl::pin<Metavariables>, tmpl::_1>>;
} // namespace control_system
