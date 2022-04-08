// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>

#include "ControlSystem/Systems/Translation.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/ControlSystem/SystemHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace control_system {
namespace {
using TranslationMap = domain::CoordinateMaps::TimeDependent::Translation<3>;

using CoordMap =
    domain::CoordinateMap<Frame::Grid, Frame::Inertial, TranslationMap>;

template <size_t DerivOrder>
void test_translation_control_system() {
  // Since we are only doing translation, turn off the
  // other control systems by passing 0 for their deriv orders
  using metavars = TestHelpers::MockMetavars<DerivOrder, 0, 0>;
  using translation_component = typename metavars::translation_component;
  using element_component = typename metavars::element_component;
  MAKE_GENERATOR(gen);

  // Global things
  domain::FunctionsOfTime::register_derived_with_charm();
  double initial_time = 0.0;
  const double initial_separation = 15.0;
  // This final time is chosen so that the damping timescales have adequate time
  // to reach the maximum damping timescale
  const double final_time = 500.0;

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
  const auto& init_trans_tuple = system_helper.init_trans_tuple();
  const std::string& translation_name = system_helper.translation_name();

  // Setup runner and all components
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavars>;
  MockRuntimeSystem runner{{std::move(domain)},
                           {std::move(initial_functions_of_time),
                            std::move(initial_measurement_timescales)}};
  ActionTesting::emplace_singleton_component_and_initialize<
      translation_component>(make_not_null(&runner), ActionTesting::NodeId{0},
                             ActionTesting::LocalCoreId{0}, init_trans_tuple);
  ActionTesting::emplace_array_component<element_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0);

  ActionTesting::set_phase(make_not_null(&runner), metavars::Phase::Testing);

  // Create coordinate map for mapping the translation to the "grid" frame
  // where the control system does its calculations
  TranslationMap translation_map{translation_name};

  CoordMap coord_map{translation_map};

  // Get the functions of time from the cache to use in the maps
  const auto& cache = ActionTesting::cache<element_component>(runner, 0);
  const auto& functions_of_time =
      Parallel::get<domain::Tags::FunctionsOfTime>(cache);

  const std::array<double, 3> velocity{{0.1, -0.2, 0.3}};

  const auto position_function = [&initial_separation,
                                  &velocity](const double time) {
    const std::array<double, 3> init_pos{{0.5 * initial_separation, 0.0, 0.0}};
    return std::pair<std::array<double, 3>, std::array<double, 3>>{
        -init_pos + velocity * time, init_pos + velocity * time};
  };

  // Run the actual control system test.
  system_helper.run_control_system_test(runner, final_time, make_not_null(&gen),
                                        position_function, coord_map);

  // Grab results
  const std::array<double, 3> grid_position_of_a =
      system_helper.grid_position_of_a();
  const std::array<double, 3> grid_position_of_b =
      system_helper.grid_position_of_b();

  // Our expected positions are just the initial positions
  const std::array<double, 3> expected_grid_position_of_a{
      {-0.5 * initial_separation, 0.0, 0.0}};
  const std::array<double, 3> expected_grid_position_of_b{
      {0.5 * initial_separation, 0.0, 0.0}};

  const auto& translation_f_of_t =
      dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<DerivOrder>&>(
          *functions_of_time.at(translation_name));

  const auto trans_and_2_derivs =
      translation_f_of_t.func_and_2_derivs(final_time);

  // The control system gets more accurate the longer you run for. This is
  // the accuracy we can achieve in this amount of time.
  Approx custom_approx = Approx::custom().epsilon(1.0e-12).scale(1.0);
  const DataVector expected_translation_2nd_deriv{3, 0.0};
  CHECK_ITERABLE_CUSTOM_APPROX(trans_and_2_derivs[2],
                               expected_translation_2nd_deriv, custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(trans_and_2_derivs[1],
                               array_to_datavector(velocity), custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(trans_and_2_derivs[0],
                               array_to_datavector(velocity * final_time),
                               custom_approx);

  CHECK_ITERABLE_CUSTOM_APPROX(expected_grid_position_of_a, grid_position_of_a,
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(expected_grid_position_of_b, grid_position_of_b,
                               custom_approx);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Systems.Translation",
                  "[ControlSystem][Unit]") {
  test_translation_control_system<2>();
  test_translation_control_system<3>();
}
}  // namespace control_system
