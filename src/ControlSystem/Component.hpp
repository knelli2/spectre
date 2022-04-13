// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <tuple>
#include <string>

#include "ControlSystem/Actions/Initialization.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/WriteData.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/Observer/Actions/RegisterSingleton.hpp"
#include "IO/Observer/Helpers.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Tags/ResourceInfo.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
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
  static_assert(tt::assert_conforms_to<
                ControlSystem, control_system::protocols::ControlSystem>);
  using chare_type = Parallel::Algorithms::Singleton;

  using metavariables = Metavariables;

  static std::string name() { return ControlSystem::name(); }

  using system = ControlSystem;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename metavariables::Phase, metavariables::Phase::Initialization,
          tmpl::list<
              Actions::SetupDataBox,
              control_system::Actions::Initialize<Metavariables, ControlSystem>,
              Initialization::Actions::RemoveOptionsAndTerminatePhase>>,
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Register,
          tmpl::list<observers::Actions::RegisterSingletonWithObserverWriter<
              control_system::Registration<ControlSystem>>>>>;

  using initialization_tags = tmpl::flatten<
      tmpl::list<Parallel::get_initialization_tags<
                     Parallel::get_initialization_actions_list<
                         phase_dependent_action_list>>,
                 Parallel::Tags::SingletonInfo<ControlComponent>>>;

  static void execute_next_phase(
      const typename Metavariables::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
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
}  // namespace control_system
