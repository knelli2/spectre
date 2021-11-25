// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/FunctionsOfTimeInitialize.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/OptionTags.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Time/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
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
  static constexpr size_t volume_dim = 1;

  using control_systems = tmpl::list<FakeControlSystem<1>, FakeControlSystem<2>,
                                     FakeControlSystem<3>>;
  using component_list =
      control_system::control_components<Metavariables, control_systems>;
};

struct MetavariablesNoControlSystems {
  using component_list = tmpl::list<>;
};

class TestCreator : public DomainCreator<1> {
 public:
  explicit TestCreator(const bool add_controlled)
      : add_controlled_(add_controlled) {}
  Domain<1> create_domain() const override { ERROR(""); }
  std::vector<std::array<size_t, 1>> initial_extents() const override {
    ERROR("");
  }
  std::vector<std::array<size_t, 1>> initial_refinement_levels()
      const override {
    ERROR("");
  }
  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override {
    const std::array<DataVector, 3> initial_values{{{-1.0}, {-2.0}, {-3.0}}};

    std::unordered_map<std::string, double> expiration_times{
        {"Uncontrolled", std::numeric_limits<double>::infinity()}};
    if (add_controlled_) {
      expiration_times["Controlled1"] = std::numeric_limits<double>::infinity();
      expiration_times["Controlled2"] = std::numeric_limits<double>::infinity();
      expiration_times["Controlled3"] = std::numeric_limits<double>::infinity();
      for (auto& [name, expr_time] : initial_expiration_times) {
        if (expiration_times.count(name) == 1) {
          expiration_times[name] = expr_time;
        }
      }
    }

    std::unordered_map<std::string,
                       std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
        result{};

    if (add_controlled_) {
      result.insert(
          {"Controlled1",
           std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
               initial_time, initial_values,
               expiration_times.at("Controlled1"))});
      result.insert(
          {"Controlled2",
           std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
               initial_time, initial_values,
               expiration_times.at("Controlled2"))});
      result.insert(
          {"Controlled3",
           std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
               initial_time, initial_values,
               expiration_times.at("Controlled3"))});
    }
    result.insert(
        {"Uncontrolled",
         std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
             initial_time, initial_values,
             expiration_times.at("Uncontrolled"))});
    return result;
  }

 private:
  bool add_controlled_{};
};

class BadCreator : public DomainCreator<1> {
 public:
  Domain<1> create_domain() const override { ERROR(""); }
  std::vector<std::array<size_t, 1>> initial_extents() const override {
    ERROR("");
  }
  std::vector<std::array<size_t, 1>> initial_refinement_levels()
      const override {
    ERROR("");
  }
  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override {
    TestCreator good_creator{true};
    auto functions_of_time =
        good_creator.functions_of_time(initial_expiration_times);

    // Mimick a domain creator that has improperly set the expiration time for
    // one of the functions of time
    const auto& function_to_replace = functions_of_time.begin()->second;
    functions_of_time[initial_expiration_times.begin()->first] =
        std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
            initial_time, function_to_replace->func_and_2_derivs(initial_time),
            std::numeric_limits<double>::infinity());

    return functions_of_time;
  }
};

template <size_t Index>
using OptionHolder = control_system::OptionHolder<FakeControlSystem<Index>>;

template <typename ControlSys>
using ControlSysInputs =
    control_system::OptionTags::ControlSystemInputs<ControlSys>;

void test_functions_of_time_tag() {
  INFO("Test FunctionsOfTimeInitialize tag");
  using fot_tag = control_system::Tags::FunctionsOfTimeInitialize;
  using Creator = tmpl::front<fot_tag::option_tags<Metavariables>>::type;

  const Creator creator = std::make_unique<TestCreator>(true);

  // Initial expiration times are set to be update_fraction *
  // min(current_timescale) where update_fraction is the argument to the
  // Controller. This value for the timescale was chosen to give an expiration
  // time between the two expiration times used above in the TestCreator
  const double timescale = 27.0;
  const TimescaleTuner tuner1({timescale}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.01,
                              0.99);
  const TimescaleTuner tuner2({0.1}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.01, 0.99);
  const Averager<2> averager(0.25, true);
  const double update_fraction = 0.3;
  const Controller<2> controller(update_fraction);

  OptionHolder<1> option_holder1(averager, controller, tuner1);
  OptionHolder<2> option_holder2(averager, controller, tuner1);
  OptionHolder<3> option_holder3(averager, controller, tuner2);

  const double initial_time_step = 1.0;
  fot_tag::type functions_of_time = fot_tag::create_from_options<Metavariables>(
      creator, initial_time, initial_time_step, option_holder1, option_holder2,
      option_holder3);

  CHECK(functions_of_time.at("Controlled1")->time_bounds()[1] ==
        initial_time + update_fraction * timescale);
  CHECK(functions_of_time.at("Controlled2")->time_bounds()[1] ==
        initial_time + update_fraction * timescale);
  CHECK(functions_of_time.at("Controlled3")->time_bounds()[1] ==
        initial_time + initial_time_step);
  CHECK(functions_of_time.at("Uncontrolled")->time_bounds()[1] ==
        std::numeric_limits<double>::infinity());

  static_assert(
      std::is_same_v<
          fot_tag::option_tags<Metavariables>,
          tmpl::list<
              domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
              ::OptionTags::InitialTime, ::OptionTags::InitialTimeStep,
              ControlSysInputs<FakeControlSystem<1>>,
              ControlSysInputs<FakeControlSystem<2>>,
              ControlSysInputs<FakeControlSystem<3>>>>);
}

SPECTRE_TEST_CASE("Unit.ControlSystem.Tags.FunctionsOfTimeInitialize",
                  "[ControlSystem][Unit]") {
  test_functions_of_time_tag();
}

// [[OutputRegex, is not controlling a function of time. Check that the
// DomainCreator you have chosen uses all of the control systems in the
// executable. The existing functions of time are]]
SPECTRE_TEST_CASE(
    "Unit.ControlSystem.Tags.FunctionsOfTimeInitialize.ExtraControlSystem",
    "[ControlSystem][Unit]") {
  ERROR_TEST();
  using fot_tag = control_system::Tags::FunctionsOfTimeInitialize;
  using Creator = tmpl::front<fot_tag::option_tags<Metavariables>>::type;

  const Creator creator = std::make_unique<TestCreator>(true);

  const TimescaleTuner tuner({1.0}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.01, 0.99);
  const Averager<2> averager(0.25, true);
  const double update_fraction = 0.3;
  const Controller<2> controller(update_fraction);

  OptionHolder<1> option_holder1(averager, controller, tuner);
  OptionHolder<2> option_holder2(averager, controller, tuner);
  OptionHolder<3> option_holder3(averager, controller, tuner);
  OptionHolder<4> option_holder4(averager, controller, tuner);

  const double initial_time_step = 1.0;
  fot_tag::type functions_of_time = fot_tag::create_from_options<Metavariables>(
      creator, initial_time, initial_time_step, option_holder1, option_holder2,
      option_holder3, option_holder4);
}

// [[OutputRegex, Is the DomainCreator you are using compatible with control
// systems?]]
SPECTRE_TEST_CASE(
    "Unit.ControlSystem.Tags.FunctionsOfTimeInitialize.ImproperExprTime",
    "[ControlSystem][Unit]") {
  ERROR_TEST();
  using fot_tag = control_system::Tags::FunctionsOfTimeInitialize;
  using Creator = tmpl::front<fot_tag::option_tags<Metavariables>>::type;

  const Creator creator = std::make_unique<BadCreator>();

  const TimescaleTuner tuner({1.0}, 10.0, 1.0e-3, 1.0e-2, 1.0e-4, 1.01, 0.99);
  const Averager<2> averager(0.25, true);
  const double update_fraction = 0.3;
  const Controller<2> controller(update_fraction);

  OptionHolder<1> option_holder1(averager, controller, tuner);
  OptionHolder<2> option_holder2(averager, controller, tuner);
  OptionHolder<3> option_holder3(averager, controller, tuner);

  const double initial_time_step = 1.0;
  fot_tag::type functions_of_time = fot_tag::create_from_options<Metavariables>(
      creator, initial_time, initial_time_step, option_holder1, option_holder2,
      option_holder3);
}
}  // namespace
