// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>

#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/ControlSystem/SystemHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system {
namespace {
void test_expansion_control_error() {
  // Since we are only doing expansion, turn off the
  // other control systems by passing 0 for their deriv orders
  constexpr size_t deriv_order = 2;
  using metavars = TestHelpers::MockMetavars<0, 0, deriv_order>;
  using element_component = typename metavars::element_component;
  MAKE_GENERATOR(gen);

  // Global things
  domain::FunctionsOfTime::register_derived_with_charm();
  double initial_time = 0.0;
  const double initial_separation = 15.0;

  // Set up the system helper.
  control_system::TestHelpers::SystemHelper<metavars> system_helper{};

  // Initialize everything within the system helper
  system_helper.setup_control_system_test(initial_time, initial_separation);

  // Get references to everything that was set up inside the system helper. The
  // domain and two functions of time are not const references because they need
  // to be moved into the runner
  auto& domain = system_helper.domain();
  auto& initial_functions_of_time = system_helper.initial_functions_of_time();
  auto& initial_measurement_timescales =
      system_helper.initial_measurement_timescales();
  const std::string& expansion_name = system_helper.expansion_name();

  // Setup runner and element component because it's the easiest way to get the
  // global cache
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavars>;
  MockRuntimeSystem runner{{std::move(domain)},
                           {std::move(initial_functions_of_time),
                            std::move(initial_measurement_timescales)}};
  ActionTesting::emplace_array_component<element_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0);
  const auto& cache = ActionTesting::cache<element_component>(runner, 0);
  const auto& functions_of_time =
      Parallel::get<domain::Tags::FunctionsOfTime>(cache);

  using QueueTuple = tuples::tagged_tuple_from_typelist<tmpl::list<
      control_system::QueueTags::Center<control_system::ah::HorizonLabel::AhA>,
      control_system::QueueTags::Center<
          control_system::ah::HorizonLabel::AhB>>>;

  // Create fake measurements. For expansion we only care about the x component
  // because that's all that is used. B is on the positive x-axis, A is on the
  // negative x-axis
  std::uniform_real_distribution<double> dist{5.0, 15.0};
  const double pos_A_x = -dist(gen);
  const double pos_B_x = dist(gen);
  QueueTuple fake_measurement_tuple{DataVector{pos_A_x, 0.0, 0.0},
                                    DataVector{pos_B_x, 0.0, 0.0}};

  using expansion_system = typename metavars::expansion_system;
  using ControlError = expansion_system::control_error;

  // This is before the first expiration time
  const double check_time = 0.1;
  const DataVector control_error =
      ControlError{}(cache, check_time, expansion_name, fake_measurement_tuple);

  const auto& expansion_f_of_t =
      dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>&>(
          *functions_of_time.at(expansion_name));
  const auto& domain_from_cache = get<domain::Tags::Domain<3>>(cache);
  // Since we haven't updated, the expansion factor should just be 1.0
  const double exp_factor = expansion_f_of_t.func(check_time)[0][0];
  const double pos_diff = pos_B_x - pos_A_x;
  const double grid_diff = domain_from_cache.excision_spheres()
                               .at("ObjectBExcisionSphere")
                               .center()[0] -
                           domain_from_cache.excision_spheres()
                               .at("ObjectAExcisionSphere")
                               .center()[0];

  const DataVector expected_control_error{exp_factor *
                                          (pos_diff / grid_diff - 1.0)};

  CHECK(control_error == expected_control_error);
}

SPECTRE_TEST_CASE("Unit.ControlSystem.ControlErrors.Expansion",
                  "[ControlSystem][Unit]") {
  test_expansion_control_error();
}
}  // namespace
}  // namespace control_system
