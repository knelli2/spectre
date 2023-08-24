// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <tuple>

#include "ControlSystem/Actions/Initialization.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \ingroup ControlSystemGroup
 * \brief The singleton parallel component responsible for managing a single
 * control system.
 *
 * The control system that this component is responsible for is specified by the
 * `ControlSystem` template parameter. This control system must conform to the
 * `control_system::Protocols::ControlSystem` protocol.
 */
template <class Metavariables, typename ControlSystem>
struct ControlComponent {
  static_assert(tt::assert_conforms_to_v<
                ControlSystem, control_system::protocols::ControlSystem>);
  using chare_type = Parallel::Algorithms::Singleton;
  using control_system = ControlSystem;

  static std::string name() { return ControlSystem::name(); }

  using metavariables = Metavariables;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<Initialization::Actions::InitializeItems<
                     ::control_system::Actions::Initialize<Metavariables,
                                                           ControlSystem>>,
                 Parallel::Actions::TerminatePhase>>>;

  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<
        ControlComponent<Metavariables, ControlSystem>>(local_cache)
        .start_phase(next_phase);
  }
};

namespace control_system {
/// \ingroup ControlSystemGroup
/// List of control components to be added to the component list of the
/// metavars
template <typename Metavariables, typename ControlSystems>
using control_components = tmpl::transform<
    ControlSystems,
    tmpl::bind<ControlComponent, tmpl::pin<Metavariables>, tmpl::_1>>;
}  // namespace control_system
