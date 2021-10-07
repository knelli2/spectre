// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>

#include "ControlSystem/DataVectorHelpers.hpp"
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
void test_rotation_control_error() {
  // Since we are only doing rotation, turn off the
  // other control systems by passing 0 for their deriv orders
  constexpr size_t deriv_order = 2;
  using metavars = TestHelpers::MockMetavars<0, deriv_order, 0>;
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
  const std::string& rotation_name = system_helper.rotation_name();

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

  using QueueTuple = tuples::tagged_tuple_from_typelist<tmpl::list<
      control_system::QueueTags::Center<control_system::ah::HorizonLabel::AhA>,
      control_system::QueueTags::Center<
          control_system::ah::HorizonLabel::AhB>>>;

  // Create fake measurements.
  std::uniform_real_distribution<double> x_dist{5.0, 15.0};
  std::uniform_real_distribution<double> yz_dist{-0.5, 0.5};
  const DataVector pos_A{{-x_dist(gen), yz_dist(gen), yz_dist(gen)}};
  const DataVector pos_B{{x_dist(gen), yz_dist(gen), yz_dist(gen)}};
  QueueTuple fake_measurement_tuple{pos_A, pos_B};

  using rotation_system = typename metavars::rotation_system;
  using ControlError = rotation_system::control_error;

  // This is before the first expiration time
  const double check_time = 0.1;
  const DataVector control_error =
      ControlError{}(cache, check_time, rotation_name, fake_measurement_tuple);

  const auto& domain_from_cache = get<domain::Tags::Domain<3>>(cache);
  const DataVector pos_diff = pos_B - pos_A;
  const DataVector grid_diff =
      array_to_datavector(domain_from_cache.excision_spheres()
                              .at("ObjectBExcisionSphere")
                              .center()) -
      array_to_datavector(domain_from_cache.excision_spheres()
                              .at("ObjectAExcisionSphere")
                              .center());

  const double grid_dot_pos = dot(grid_diff, pos_diff);
  const DataVector grid_cross_pos = cross(grid_diff, pos_diff);

  const DataVector expected_control_error = grid_cross_pos / grid_dot_pos;

  CHECK(control_error == expected_control_error);
}

SPECTRE_TEST_CASE("Unit.ControlSystem.ControlErrors.Rotation",
                  "[ControlSystem][Unit]") {
  test_rotation_control_error();
}
}  // namespace
}  // namespace control_system
