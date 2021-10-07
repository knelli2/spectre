// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>

#include "ControlSystem/ControlErrors/Expansion.hpp"
#include "ControlSystem/ControlErrors/Rotation.hpp"
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
void test_translation_control_error() {
  // Since we are only doing translation, turn off the
  // other control systems by passing 0 for their deriv orders
  constexpr size_t deriv_order = 2;
  using metavars = TestHelpers::MockMetavars<deriv_order, 0, 0>;
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
  const std::string& translation_name = system_helper.translation_name();
  const std::string& rotation_name = system_helper.rotation_name();
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

  // Create fake measurements.
  std::uniform_real_distribution<double> x_dist{5.0, 15.0};
  std::uniform_real_distribution<double> yz_dist{-0.5, 0.5};
  const DataVector pos_A{{-x_dist(gen), yz_dist(gen), yz_dist(gen)}};
  const DataVector pos_B{{x_dist(gen), yz_dist(gen), yz_dist(gen)}};
  QueueTuple fake_measurement_tuple{pos_A, pos_B};

  using translation_system = typename metavars::translation_system;
  using ControlError = translation_system::control_error;

  // This is before the first expiration time
  const double check_time = 0.1;
  const DataVector control_error = ControlError{}(
      cache, check_time, translation_name, fake_measurement_tuple);

  const DataVector rotation_control_error =
      control_system::ControlErrors::Rotation{}(
          cache, check_time, rotation_name, fake_measurement_tuple);
  const double expansion_control_error =
      control_system::ControlErrors::Expansion{}(
          cache, check_time, expansion_name, fake_measurement_tuple)[0];

  // deriv_order is 0 here because we don't need a real expansion function of
  // time for translation. Inside SystemHelper, a PiecewisePolynomial<0> with
  // the expansion factor = 1 is added
  const auto& expansion_f_of_t =
      dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<0>&>(
          *functions_of_time.at(expansion_name));
  const double exp_factor = expansion_f_of_t.func(check_time)[0][0];
  const auto& domain_from_cache = get<domain::Tags::Domain<3>>(cache);
  const DataVector grid_B =
      array_to_datavector(domain_from_cache.excision_spheres()
                              .at("ObjectBExcisionSphere")
                              .center());

  const DataVector rot_control_err_cross_grid =
      cross(rotation_control_error, grid_B);

  // The quaternion should be the unit quaternion (1,0,0,0) which means the
  // quaternion multiplication in the translation control error is the identity
  // so we avoid actually doing quaternion multiplication
  const DataVector expected_control_error =
      exp_factor * (pos_B - grid_B - rot_control_err_cross_grid -
                    expansion_control_error / exp_factor * grid_B);

  Approx custom_approx = Approx::custom().epsilon(1.0e-14).scale(1.0);
  CHECK_ITERABLE_CUSTOM_APPROX(control_error, expected_control_error,
                               custom_approx);
}

SPECTRE_TEST_CASE("Unit.ControlSystem.ControlErrors.Translation",
                  "[ControlSystem][Unit]") {
  test_translation_control_error();
}
}  // namespace
}  // namespace control_system
