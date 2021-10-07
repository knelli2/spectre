// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Expansion.hpp"
#include "ControlSystem/ControlErrors/Rotation.hpp"
#include "ControlSystem/DataVectorHelpers.hpp"
#include "ControlSystem/Systems/Expansion.hpp"
#include "ControlSystem/Systems/Rotation.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestingFramework.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Utilities/TMPL.hpp"

namespace control_system::TestHelpers {
template <size_t DerivOrder>
struct MockExpansionSystem : control_system::Expansion<DerivOrder> {
  struct MeasurementQueue : db::SimpleTag {
    using type = LinkedMessageQueue<
        double,
        tmpl::list<control_system::QueueTags::Center<ah::HorizonLabel::AhA>,
                   control_system::QueueTags::Center<ah::HorizonLabel::AhB>>>;
  };

  static constexpr size_t deriv_order = DerivOrder;

  struct process_measurement {
    template <typename Component, ah::HorizonLabel Horizon,
              typename Metavariables>
    static void apply(
        gsl::not_null<ActionTesting::MockRuntimeSystem<Metavariables>*> runner,
        const std::array<double, 3>& center_array,
        const LinkedMessageId<double>& measurement_id) {
      const DataVector center = array_to_datavector(center_array);

      ActionTesting::simple_action<
          Component,
          ::Actions::UpdateMessageQueue<
              control_system::QueueTags::Center<Horizon>, MeasurementQueue,
              control_system::UpdateControlSystem<
                  deriv_order, control_system::ControlErrors::Expansion>>>(
          runner, 0, measurement_id, center);
    }
  };
};

template <size_t DerivOrder>
struct MockRotationSystem : control_system::Rotation<DerivOrder> {
  struct MeasurementQueue : db::SimpleTag {
    using type = LinkedMessageQueue<
        double,
        tmpl::list<control_system::QueueTags::Center<ah::HorizonLabel::AhA>,
                   control_system::QueueTags::Center<ah::HorizonLabel::AhB>>>;
  };

  static constexpr size_t deriv_order = DerivOrder;

  struct process_measurement {
    template <typename Component, ah::HorizonLabel Horizon,
              typename Metavariables>
    static void apply(
        gsl::not_null<ActionTesting::MockRuntimeSystem<Metavariables>*> runner,
        const std::array<double, 3>& center_array,
        const LinkedMessageId<double>& measurement_id) {
      const DataVector center = array_to_datavector(center_array);

      ActionTesting::simple_action<
          Component,
          ::Actions::UpdateMessageQueue<
              control_system::QueueTags::Center<Horizon>, MeasurementQueue,
              control_system::UpdateControlSystem<
                  deriv_order, control_system::ControlErrors::Rotation>>>(
          runner, 0, measurement_id, center);
    }
  };
};

template <typename ControlSystem>
using init_simple_tags =
    tmpl::list<control_system::Tags::Averager<ControlSystem::deriv_order - 1>,
               control_system::Tags::TimescaleTuner,
               control_system::Tags::ControlSystemName,
               control_system::Tags::Controller<ControlSystem::deriv_order>,
               typename ControlSystem::MeasurementQueue>;

template <typename Metavariables, typename MockControlSystem>
struct MockControlComponent {
  using array_index = int;
  using component_being_mocked =
      ControlComponent<Metavariables, MockControlSystem>;
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockSingletonChare;

  using simple_tags = init_simple_tags<MockControlSystem>;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename metavariables::Phase, metavariables::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>>;
};

template <typename Metavars>
struct MockElementComponent {
  using array_index = int;
  using chare_type = ActionTesting::MockArrayChare;

  using metavariables = Metavars;

  using simple_tags = tmpl::list<>;

  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Metavars::volume_dim>>;

  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize,
                 control_system::Tags::MeasurementTimescales>;

  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename metavariables::Phase,
                                        metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};
}  // namespace control_system::TestHelpers
