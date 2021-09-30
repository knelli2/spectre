// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <memory>
#include <utility>

#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Trigger.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/MockDistributedObject.hpp"
#include "Framework/MockRuntimeSystem.hpp"
#include "Framework/MockRuntimeSystemFreeFunctions.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct LabelA {};
struct LabelB {};

using measurement = control_system::TestHelpers::Measurement<LabelA>;
using SystemA = control_system::TestHelpers::System<LabelA, measurement>;
using SystemB = control_system::TestHelpers::System<LabelB, measurement>;

using MeasureTrigger = control_system::Trigger<tmpl::list<SystemA, SystemB>>;

using MeasurementFoT = domain::FunctionsOfTime::PiecewisePolynomial<0>;

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;

  using simple_tags = tmpl::list<Tags::Time>;
  using compute_tags = tmpl::list<Parallel::Tags::FromGlobalCache<
      control_system::Tags::MeasurementTimescales>>;
  using mutable_global_cache_tags =
      tmpl::list<control_system::Tags::MeasurementTimescales>;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags, compute_tags>>>>;
};

struct Metavariables {
  using component_list = tmpl::list<Component<Metavariables>>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<DenseTrigger, tmpl::list<MeasureTrigger>>>;
  };

  enum class Phase { Initialization, Testing, Exit };
};
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Trigger", "[Domain][Unit]") {
  Parallel::register_classes_with_charm<MeasurementFoT>();
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<Metavariables>;
  using component = Component<Metavariables>;
  const component* const component_p = nullptr;

  control_system::Tags::MeasurementTimescales::type measurement_timescales{};
  measurement_timescales["LabelA"] =
      std::unique_ptr<MeasurementFoT>{new MeasurementFoT(0.0, {{{3.0}}}, 1.5)};
  measurement_timescales["LabelB"] =
      std::unique_ptr<MeasurementFoT>{new MeasurementFoT(0.0, {{{1.0}}}, 0.5)};
  measurement_timescales["DifferentMeasurement"] =
      std::unique_ptr<MeasurementFoT>{new MeasurementFoT(0.0, {{{1.0}}}, 0.1)};

  MockRuntimeSystem runner{{}, {std::move(measurement_timescales)}};
  ActionTesting::emplace_array_component_and_initialize<component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0, {0.0});
  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::Testing);

  auto& box = ActionTesting::get_databox<
      component,
      tmpl::list<Tags::Time, Parallel::Tags::FromGlobalCache<
                                 control_system::Tags::MeasurementTimescales>>>(
      make_not_null(&runner), 0);
  auto& cache = ActionTesting::cache<component>(runner, 0);

  const auto set_time = [&box](const double time) noexcept {
    db::mutate<Tags::Time>(
        make_not_null(&box),
        [&time](const gsl::not_null<double*> box_time) noexcept {
          *box_time = time;
        });
  };

  MeasureTrigger typed_trigger{};
  DenseTrigger& trigger = typed_trigger;

  CHECK(trigger.is_ready(box, cache, 0, component_p));
  {
    const auto is_triggered = trigger.is_triggered(box);
    CHECK(is_triggered.is_triggered);
    CHECK(is_triggered.next_check == 1.0);
  }

  set_time(0.25);

  CHECK(trigger.is_ready(box, cache, 0, component_p));
  {
    const auto is_triggered = trigger.is_triggered(box);
    CHECK(not is_triggered.is_triggered);
    CHECK(is_triggered.next_check == 1.0);
  }

  set_time(0.75);

  CHECK(not trigger.is_ready(box, cache, 0, component_p));

  struct UpdateMeasurementTimescale {
    static void apply(
        const gsl::not_null<control_system::Tags::MeasurementTimescales::type*>
            timescales) noexcept {
      static_cast<MeasurementFoT&>(*timescales->at("LabelB"))
          .update(0.5, {4.0}, 4.0);
    }
  };

  Parallel::mutate<control_system::Tags::MeasurementTimescales,
                   UpdateMeasurementTimescale>(cache);

  CHECK(trigger.is_ready(box, cache, 0, component_p));
  {
    const auto is_triggered = trigger.is_triggered(box);
    CHECK(not is_triggered.is_triggered);
    CHECK(is_triggered.next_check == 1.0);
  }

  set_time(1.0);

  CHECK(trigger.is_ready(box, cache, 0, component_p));
  {
    const auto is_triggered = trigger.is_triggered(box);
    CHECK(is_triggered.is_triggered);
    CHECK(is_triggered.next_check == 4.0);
  }
}
