// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <string>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataVector.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
const double initial_time = 2.0;

template <size_t Index>
struct FakeControlSystem
    : tt::ConformsTo<control_system::protocols::ControlSystem> {
  static constexpr size_t deriv_order = 2;
  static std::string name() { return "Controlled"s + get_output(Index); }
  using measurement = control_system::TestHelpers::Measurement<
      control_system::TestHelpers::TestStructs_detail::LabelA>;
  using simple_tags = tmpl::list<>;
  struct process_measurement {
    using argument_tags = tmpl::list<>;
  };
};

struct Metavariables {
  using control_systems = tmpl::list<FakeControlSystem<1>, FakeControlSystem<2>,
                                     FakeControlSystem<3>>;
  using component_list =
      control_system::control_components<Metavariables, control_systems>;
};

struct MetavariablesNoControlSystems {
  using component_list = tmpl::list<>;
};

template <size_t Index>
using OptionHolder = control_system::OptionHolder<FakeControlSystem<Index>>;

void test_measurement_tag() {
  INFO("Test measurement tag");
  using measurement_tag = control_system::Tags::MeasurementTimescales;
  static_assert(
      tmpl::size<
          measurement_tag::option_tags<MetavariablesNoControlSystems>>::value ==
      2);
  static_assert(
      tmpl::size<measurement_tag::option_tags<Metavariables>>::value == 5);

  const double time_step = 0.2;
  {
    const double timescale1 = 27.0;
    const TimescaleTuner tuner1({timescale1}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4,
                                1.01, 0.99);
    const double timescale2 = 0.5;
    const TimescaleTuner tuner2({timescale2}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4,
                                1.01, 0.99);
    const double averaging_fraction = 0.25;
    const Averager<2> averager(averaging_fraction, true);
    const double update_fraction = 0.3;
    const Controller<2> controller(update_fraction);

    OptionHolder<1> option_holder1(averager, controller, tuner1);
    OptionHolder<2> option_holder2(averager, controller, tuner1);
    OptionHolder<3> option_holder3(averager, controller, tuner2);

    const measurement_tag::type timescales =
        measurement_tag::create_from_options<Metavariables>(
            initial_time, time_step, option_holder1, option_holder2,
            option_holder3);
    CHECK(timescales.size() == 3);
    // The lack of expiration is a placeholder until the control systems
    // have been implemented sufficiently to manage their timescales.
    const double expr_time1 = initial_time + update_fraction * timescale1;
    const double expr_time2 = initial_time + time_step;
    const double measure_time1 =
        averaging_fraction * (expr_time1 - initial_time);
    const double measure_time2 = time_step;
    CHECK(timescales.at("Controlled1")->time_bounds() ==
          std::array{initial_time, expr_time1});
    CHECK(timescales.at("Controlled1")->func(2.0)[0] ==
          DataVector{measure_time1});
    CHECK(timescales.at("Controlled1")->func(3.0)[0] ==
          DataVector{measure_time1});
    CHECK(timescales.at("Controlled2")->time_bounds() ==
          std::array{initial_time, expr_time1});
    CHECK(timescales.at("Controlled2")->func(2.0)[0] ==
          DataVector{measure_time1});
    CHECK(timescales.at("Controlled2")->func(3.0)[0] ==
          DataVector{measure_time1});
    CHECK(timescales.at("Controlled3")->time_bounds() ==
          std::array{initial_time, expr_time2});
    CHECK(timescales.at("Controlled3")->func(2.1)[0] ==
          DataVector{measure_time2});
  }
  {
    // Verify that no control systems means no measurement timescales
    const measurement_tag::type timescales =
        measurement_tag::create_from_options<MetavariablesNoControlSystems>(
            initial_time, time_step);
    CHECK(timescales.empty());
  }
  {
    // Verify that negative time steps are accepted with no control
    // systems.
    const measurement_tag::type timescales =
        measurement_tag::create_from_options<MetavariablesNoControlSystems>(
            initial_time, -time_step);
    CHECK(timescales.empty());
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Tags.MeasurementTimescales",
                  "[ControlSystem][Unit]") {
  test_measurement_tag();
}

// [[OutputRegex, Control systems can only be used in forward-in-time
// evolutions.]]
SPECTRE_TEST_CASE("Unit.ControlSystem.Tags.MeasurementTimescales.Backwards",
                  "[ControlSystem][Unit]") {
  ERROR_TEST();
  const TimescaleTuner tuner1({27.0}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.01, 0.99);
  const TimescaleTuner tuner2({0.1}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.01, 0.99);
  const Averager<2> averager(0.25, true);
  const Controller<2> controller(0.3);

  OptionHolder<1> option_holder1(averager, controller, tuner1);
  OptionHolder<2> option_holder2(averager, controller, tuner1);
  OptionHolder<3> option_holder3(averager, controller, tuner2);

  using measurement_tag = control_system::Tags::MeasurementTimescales;
  measurement_tag::create_from_options<Metavariables>(
      initial_time, -1.0, option_holder1, option_holder2, option_holder3);
}
