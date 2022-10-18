// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <limits>
#include <memory>
#include <utility>

#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/Trigger.hpp"
#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/MockDistributedObject.hpp"
#include "Framework/MockRuntimeSystem.hpp"
#include "Framework/MockRuntimeSystemFreeFunctions.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct LabelA {};
struct LabelB {};
struct LabelC {};

using measurement = control_system::TestHelpers::Measurement<LabelA>;
using SystemA = control_system::TestHelpers::System<2, LabelA, measurement>;
using SystemB = control_system::TestHelpers::System<2, LabelB, measurement>;
using SystemC = control_system::TestHelpers::System<2, LabelC, measurement>;

using MeasureTrigger =
    control_system::Trigger<tmpl::list<SystemA, SystemB, SystemC>>;

using MeasurementFoT = domain::FunctionsOfTime::PiecewisePolynomial<0>;

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;

  using simple_tags = tmpl::list<Tags::Time>;
  using compute_tags =
      tmpl::list<Parallel::Tags::FromGlobalCache<
                     ::domain::Tags::FunctionsOfTimeInitialize>,
                 Parallel::Tags::FromGlobalCache<
                     control_system::Tags::MeasurementTimescales>>;
  using mutable_global_cache_tags =
      tmpl::list<::domain::Tags::FunctionsOfTimeInitialize,
                 control_system::Tags::MeasurementTimescales>;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags, compute_tags>>>>;
};

struct Metavariables {
  using component_list = tmpl::list<Component<Metavariables>>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<DenseTrigger, tmpl::list<MeasureTrigger>>>;
  };
};

using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<Metavariables>;
using component = Component<Metavariables>;

void test_trigger_no_replace() {
  Parallel::register_classes_with_charm<MeasurementFoT>();
  const component* const component_p = nullptr;

  control_system::Tags::MeasurementTimescales::type measurement_timescales{};
  measurement_timescales["LabelA"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{3.0}}, 1.5);
  measurement_timescales["LabelB"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{1.0}}, 0.6);
  measurement_timescales["LabelC"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{0.5}}, 0.5);
  measurement_timescales["DifferentMeasurement"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{1.0}}, 0.1);
  ::domain::Tags::FunctionsOfTimeInitialize::type functions_of_time_init{};
  const double one_plus_eps = 1.0 + std::numeric_limits<double>::epsilon();
  functions_of_time_init["FakeFoT"] = std::make_unique<MeasurementFoT>(
      0.0, std::array{DataVector{1.0}}, one_plus_eps);

  MockRuntimeSystem runner{
      {},
      {std::move(functions_of_time_init), std::move(measurement_timescales)}};
  ActionTesting::emplace_array_component_and_initialize<component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0, {0.0});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto& box = ActionTesting::get_databox<component>(make_not_null(&runner), 0);
  auto& cache = ActionTesting::cache<component>(runner, 0);

  const auto set_time = [&box](const double time) {
    db::mutate<Tags::Time>(
        make_not_null(&box),
        [&time](const gsl::not_null<double*> box_time) { *box_time = time; });
  };

  MeasureTrigger typed_trigger = serialize_and_deserialize(MeasureTrigger{});
  DenseTrigger& trigger = typed_trigger;

  // At the initial time, the trigger should be triggered and we should know
  // the next check time
  CHECK(trigger.is_triggered(box, cache, 0, component_p) ==
        std::optional{true});
  CHECK(trigger.next_check_time(box, cache, 0, component_p) ==
        std::optional{0.5});

  set_time(0.25);

  // Set the time to sometime before the next check time. Shouldn't be
  // triggered, but should still have the same check time as before
  CHECK(trigger.is_triggered(box, cache, 0, component_p) ==
        std::optional{false});
  CHECK(trigger.next_check_time(box, cache, 0, component_p) ==
        std::optional{0.5});

  set_time(0.5);

  // Now at the previous next check time, we should trigger. The new next check
  // time should be within roundoff of the minimum function of time expiration
  // time, so the next check time will be set to exactly the minimum expiration
  // time
  CHECK(trigger.is_triggered(box, cache, 0, component_p) ==
        std::optional{true});
  CHECK(trigger.next_check_time(box, cache, 0, component_p) ==
        std::optional{one_plus_eps});

  set_time(0.75);

  // Another intermediate time where we shouldn't trigger. At this point, some
  // of the measurement timescales have expired and have not been updated yet,
  // so we cannot calculate the next check time. It should be nullopt
  CHECK(trigger.is_triggered(box, cache, 0, component_p) ==
        std::optional{false});
  CHECK(not trigger.next_check_time(box, cache, 0, component_p).has_value());

  // Update the measurement timescales
  Parallel::mutate<control_system::Tags::MeasurementTimescales,
                   control_system::UpdateFunctionOfTime>(cache, "LabelB"s, 0.6,
                                                         DataVector{2.0}, 4.0);
  Parallel::mutate<control_system::Tags::MeasurementTimescales,
                   control_system::UpdateFunctionOfTime>(cache, "LabelC"s, 0.5,
                                                         DataVector{1.0}, 4.0);

  // Now we should be able to calculate the next check time once again and it
  // should be the same as it was before, since the current time hasn't changed.
  CHECK(trigger.next_check_time(box, cache, 0, component_p) ==
        std::optional{one_plus_eps});

  set_time(one_plus_eps);

  // Now the time is at the minimum function of time expiration time, so we
  // trigger. However, the functions of time haven't been updated yet so we
  // can't determine the next check time
  CHECK(trigger.is_triggered(box, cache, 0, component_p) ==
        std::optional{true});
  CHECK(not trigger.next_check_time(box, cache, 0, component_p).has_value());

  // Update the functions of time
  Parallel::mutate<domain::Tags::FunctionsOfTime,
                   control_system::UpdateFunctionOfTime>(
      cache, "FakeFoT"s, one_plus_eps, DataVector{2.0}, 4.0);

  // Now we can calculate the next check time
  CHECK(trigger.next_check_time(box, cache, 0, component_p) ==
        std::optional{2.0});

  // Update a measurement timescale that is now expired
  Parallel::mutate<control_system::Tags::MeasurementTimescales,
                   control_system::UpdateFunctionOfTime>(cache, "LabelA"s, 1.5,
                                                         DataVector{2.0}, 4.0);

  set_time(2.0);

  // At the next trigger time and all the expiration times are ok so we can
  // calculate the next check time
  CHECK(trigger.is_triggered(box, cache, 0, component_p) ==
        std::optional{true});
  CHECK(trigger.next_check_time(box, cache, 0, component_p) ==
        std::optional{3.0});
}

void test_trigger_with_replace() {
  Parallel::register_classes_with_charm<MeasurementFoT>();
  const component* const component_p = nullptr;

  control_system::Tags::MeasurementTimescales::type measurement_timescales{};
  measurement_timescales["LabelA"] = std::make_unique<MeasurementFoT>(
      0.0, std::array{DataVector{std::numeric_limits<double>::infinity()}},
      std::numeric_limits<double>::infinity());
  measurement_timescales["LabelB"] = std::make_unique<MeasurementFoT>(
      0.0, std::array{DataVector{std::numeric_limits<double>::infinity()}},
      std::numeric_limits<double>::infinity());
  measurement_timescales["LabelC"] = std::make_unique<MeasurementFoT>(
      0.0, std::array{DataVector{std::numeric_limits<double>::infinity()}},
      std::numeric_limits<double>::infinity());
  ::domain::Tags::FunctionsOfTimeInitialize::type functions_of_time_init{};
  functions_of_time_init["FakeFoT"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{1.0}},
                                       std::numeric_limits<double>::infinity());

  MockRuntimeSystem runner{
      {},
      {std::move(functions_of_time_init), std::move(measurement_timescales)}};
  ActionTesting::emplace_array_component_and_initialize<component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0, {0.0});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto& box = ActionTesting::get_databox<component>(make_not_null(&runner), 0);
  auto& cache = ActionTesting::cache<component>(runner, 0);

  MeasureTrigger typed_trigger = serialize_and_deserialize(MeasureTrigger{});
  DenseTrigger& trigger = typed_trigger;

  const auto is_triggered = trigger.is_triggered(box, cache, 0, component_p);
  CHECK(is_triggered == std::optional{false});
  const auto next_check = trigger.next_check_time(box, cache, 0, component_p);
  CHECK(next_check == std::optional{std::numeric_limits<double>::infinity()});
}

void test_fot_expr_infinity() {
  Parallel::register_classes_with_charm<MeasurementFoT>();
  const component* const component_p = nullptr;

  control_system::Tags::MeasurementTimescales::type measurement_timescales{};
  measurement_timescales["LabelA"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{0.1}},
                                       std::numeric_limits<double>::infinity());
  measurement_timescales["LabelB"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{0.2}},
                                       std::numeric_limits<double>::infinity());
  measurement_timescales["LabelC"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{0.3}},
                                       std::numeric_limits<double>::infinity());
  ::domain::Tags::FunctionsOfTimeInitialize::type functions_of_time_init{};
  functions_of_time_init["FakeFoT"] =
      std::make_unique<MeasurementFoT>(0.0, std::array{DataVector{1.0}},
                                       std::numeric_limits<double>::infinity());

  MockRuntimeSystem runner{
      {},
      {std::move(functions_of_time_init), std::move(measurement_timescales)}};
  ActionTesting::emplace_array_component_and_initialize<component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0, {0.0});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto& box = ActionTesting::get_databox<component>(make_not_null(&runner), 0);
  auto& cache = ActionTesting::cache<component>(runner, 0);

  MeasureTrigger typed_trigger = serialize_and_deserialize(MeasureTrigger{});
  DenseTrigger& trigger = typed_trigger;

  const auto is_triggered = trigger.is_triggered(box, cache, 0, component_p);
  CHECK(is_triggered == std::optional{true});
  const auto next_check = trigger.next_check_time(box, cache, 0, component_p);
  CHECK(next_check == std::optional{0.1});
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Trigger", "[Domain][Unit]") {
  test_trigger_no_replace();
  test_trigger_with_replace();
  test_fot_expr_infinity();
}
