// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/FunctionsOfTimeInitialize.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "Framework/TestCreation.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"

namespace {
struct LabelA {};
using system = control_system::TestHelpers::System<
    2, LabelA, control_system::TestHelpers::Measurement<LabelA>>;

struct Metavars {
  static constexpr bool override_functions_of_time = true;
};
struct MetavarsEmpty {};

void test_all_tags() {
  INFO("Test all tags");
  using write_tag = control_system::Tags::WriteDataToDisk;
  TestHelpers::db::test_simple_tag<write_tag>("WriteDataToDisk");
  using averager_tag = control_system::Tags::Averager<system>;
  TestHelpers::db::test_simple_tag<averager_tag>("Averager");
  using timescaletuner_tag = control_system::Tags::TimescaleTuner<system>;
  TestHelpers::db::test_simple_tag<timescaletuner_tag>("TimescaleTuner");
  using controller_tag = control_system::Tags::Controller<system>;
  TestHelpers::db::test_simple_tag<controller_tag>("Controller");
  using fot_tag = control_system::Tags::FunctionsOfTimeInitialize;
  TestHelpers::db::test_simple_tag<fot_tag>("FunctionsOfTime");

  using control_error_tag = control_system::Tags::ControlError<system>;
  TestHelpers::db::test_simple_tag<control_error_tag>("ControlError");

  using active_tag = control_system::Tags::IsActive<system>;
  TestHelpers::db::test_simple_tag<active_tag>("IsActive");
  // For some reason, clang-10 says this const static variable is unused so
  // just static assert it here.
  static_assert(Metavars::override_functions_of_time);
  CHECK(active_tag::create_from_options<MetavarsEmpty>());
  CHECK_FALSE(active_tag::create_from_options<Metavars>(
      {"/fake/path"}, {{"FakeSpecName", "LabelA"}}));
  CHECK(active_tag::create_from_options<Metavars>({"/fake/path"}, {}));
  CHECK(active_tag::create_from_options<Metavars>(
      std::nullopt, {{"FakeSpecName", "LabelA"}}));

  using measurement_tag = control_system::Tags::MeasurementTimescales;
  TestHelpers::db::test_simple_tag<measurement_tag>("MeasurementTimescales");
}

void test_control_sys_inputs() {
  INFO("Test control system inputs");
  const double decrease_timescale_threshold = 1.0e-2;
  const double increase_timescale_threshold = 1.0e-4;
  const double increase_factor = 1.01;
  const double decrease_factor = 0.99;
  const double max_timescale = 10.0;
  const double min_timescale = 1.0e-3;
  const TimescaleTuner expected_tuner(
      {1.}, max_timescale, min_timescale, decrease_timescale_threshold,
      increase_timescale_threshold, increase_factor, decrease_factor);
  const Averager<1> expected_averager(0.25, true);
  const Controller<2> expected_controller(0.3);
  const std::string expected_name{"LabelA"};

  const auto input_holder = TestHelpers::test_option_tag<
      control_system::OptionTags::ControlSystemInputs<system>>(
      "Averager:\n"
      "  AverageTimescaleFraction: 0.25\n"
      "  Average0thDeriv: true\n"
      "Controller:\n"
      "  UpdateFraction: 0.3\n"
      "TimescaleTuner:\n"
      "  InitialTimescales: [1.]\n"
      "  MinTimescale: 1e-3\n"
      "  MaxTimescale: 10.\n"
      "  DecreaseThreshold: 1e-2\n"
      "  IncreaseThreshold: 1e-4\n"
      "  IncreaseFactor: 1.01\n"
      "  DecreaseFactor: 0.99\n"
      "ControlError:\n");
  CHECK(expected_averager == input_holder.averager);
  CHECK(expected_controller == input_holder.controller);
  CHECK(expected_tuner == input_holder.tuner);
  CHECK(expected_name ==
        std::decay_t<decltype(input_holder)>::control_system::name());

  const auto write_data =
      TestHelpers::test_option_tag<control_system::OptionTags::WriteDataToDisk>(
          "true");
  CHECK(write_data);
  // We don't check the control error because the example one is empty and
  // doesn't have a comparison operator. Once a control error is added that
  // contains member data (and thus, options), then it can be tested
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Tags", "[ControlSystem][Unit]") {
  test_all_tags();
  test_control_sys_inputs();
}
