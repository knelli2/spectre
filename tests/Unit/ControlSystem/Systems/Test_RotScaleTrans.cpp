// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <string>

#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/ControlSystem/SystemHelpers.hpp"
#include "Helpers/PointwiseFunctions/PostNewtonian/BinaryTrajectories.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace control_system {
namespace {
using TranslationMap = domain::CoordinateMaps::TimeDependent::Translation<3>;
using RotationMap = domain::CoordinateMaps::TimeDependent::Rotation<3>;
using ExpansionMap = domain::CoordinateMaps::TimeDependent::CubicScale<3>;

using CoordMap =
    domain::CoordinateMap<Frame::Grid, Frame::Inertial, ExpansionMap,
                          RotationMap, TranslationMap>;

template <size_t TranslationDerivOrder, size_t RotationDerivOrder,
          size_t ExpansionDerivOrder>
void test_rotscaletrans_control_system(const double rotation_eps = 5.0e-5) {
  using metavars =
      TestHelpers::MockMetavars<TranslationDerivOrder, RotationDerivOrder,
                                ExpansionDerivOrder>;
  using element_component = typename metavars::element_component;
  using translation_component = typename metavars::translation_component;
  using rotation_component = typename metavars::rotation_component;
  using expansion_component = typename metavars::expansion_component;
  MAKE_GENERATOR(gen);

  // Global things
  domain::FunctionsOfTime::register_derived_with_charm();
  double initial_time = 0.0;
  const double initial_separation = 15.0;
  // This final time is chosen so that the damping timescales have adequate time
  // to reach the maximum damping timescale
  const double final_time = 500.0;

  // Set up the system helper
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
  const auto& init_rot_tuple = system_helper.init_rot_tuple();
  const auto& init_exp_tuple = system_helper.init_exp_tuple();

  // Setup runner and all components
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavars>;
  MockRuntimeSystem runner{{std::move(domain)},
                           {std::move(initial_functions_of_time),
                            std::move(initial_measurement_timescales)}};
  if constexpr (TranslationDerivOrder > 0) {
    ActionTesting::emplace_singleton_component_and_initialize<
        translation_component>(make_not_null(&runner), ActionTesting::NodeId{0},
                               ActionTesting::LocalCoreId{0}, init_trans_tuple);
  }
  if constexpr (RotationDerivOrder > 0) {
    ActionTesting::emplace_singleton_component_and_initialize<
        rotation_component>(make_not_null(&runner), ActionTesting::NodeId{0},
                            ActionTesting::LocalCoreId{0}, init_rot_tuple);
  }
  if constexpr (ExpansionDerivOrder > 0) {
    ActionTesting::emplace_singleton_component_and_initialize<
        expansion_component>(make_not_null(&runner), ActionTesting::NodeId{0},
                             ActionTesting::LocalCoreId{0}, init_exp_tuple);
  }
  ActionTesting::emplace_array_component<element_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0);

  ActionTesting::set_phase(make_not_null(&runner), metavars::Phase::Testing);

  // PN orbits
  const std::array<double, 3> initial_velocity{0.1, -0.2, 0.3};
  const BinaryTrajectories binary_trajectories{initial_separation,
                                               initial_velocity};

  const std::string& translation_name = system_helper.translation_name();
  const std::string& rotation_name = system_helper.rotation_name();
  const std::string& expansion_name = system_helper.expansion_name();

  // Create coordinate maps for mapping the PN trajectories to the "grid" frame
  // where the control system does its calculations
  TranslationMap translation_map{translation_name};
  RotationMap rotation_map{rotation_name};
  // The outer boudary is at 1000.0 so that we don't have to worry about it.
  ExpansionMap expansion_map{1000.0, expansion_name,
                             expansion_name + "OuterBoundary"s};

  CoordMap coord_map{expansion_map, rotation_map, translation_map};

  // Get the functions of time from the cache to use in the maps
  const auto& cache = ActionTesting::cache<element_component>(runner, 0);
  const auto& functions_of_time =
      Parallel::get<domain::Tags::FunctionsOfTime>(cache);

  const auto position_function = [&binary_trajectories](const double time) {
    return binary_trajectories.positions(time);
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

  const auto& rotation_f_of_t = dynamic_cast<
      domain::FunctionsOfTime::QuaternionFunctionOfTime<RotationDerivOrder>&>(
      *functions_of_time.at(rotation_name));
  const auto& translation_f_of_t = dynamic_cast<
      domain::FunctionsOfTime::PiecewisePolynomial<TranslationDerivOrder>&>(
      *functions_of_time.at(translation_name));
  const auto& expansion_f_of_t = dynamic_cast<
      domain::FunctionsOfTime::PiecewisePolynomial<ExpansionDerivOrder>&>(
      *functions_of_time.at(expansion_name));

  const auto omega = rotation_f_of_t.angle_func_and_deriv(final_time)[1];
  const auto trans_and_2_derivs =
      translation_f_of_t.func_and_2_derivs(final_time);
  const auto expansion_factor = expansion_f_of_t.func(final_time)[0][0];

  // The control system gets more accurate the longer you run for. This is
  // the accuracy we can achieve in this amount of time.
  Approx custom_approx1 = Approx::custom().epsilon(5.0e-5).scale(1.0);
  Approx custom_approx2 = Approx::custom().epsilon(rotation_eps).scale(1.0);

  const DataVector expected_omega{
      {0.0, 0.0, binary_trajectories.omega(final_time)}};
  const DataVector expected_translation_2nd_deriv{3, 0.0};
  CHECK_ITERABLE_CUSTOM_APPROX(expected_omega, omega, custom_approx1);
  CHECK_ITERABLE_CUSTOM_APPROX(trans_and_2_derivs[2],
                               expected_translation_2nd_deriv, custom_approx1);
  CHECK_ITERABLE_CUSTOM_APPROX(trans_and_2_derivs[1],
                               array_to_datavector(initial_velocity),
                               custom_approx1);
  CHECK_ITERABLE_CUSTOM_APPROX(
      trans_and_2_derivs[0], array_to_datavector(initial_velocity * final_time),
      custom_approx1);
  CHECK_ITERABLE_CUSTOM_APPROX(
      binary_trajectories.separation(final_time) / initial_separation,
      expansion_factor, custom_approx1);

  CHECK_ITERABLE_CUSTOM_APPROX(expected_grid_position_of_a, grid_position_of_a,
                               custom_approx2);
  CHECK_ITERABLE_CUSTOM_APPROX(expected_grid_position_of_b, grid_position_of_b,
                               custom_approx2);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Systems.RotScaleTrans",
                  "[ControlSystem][Unit]") {
  // For PN deriv order 2, we need a different epsilon because controlling the
  // first derivative of omega (second derivative of the angle) doesn't give as
  // accurate results for the final positions of the objects in the grid frame.
  test_rotscaletrans_control_system<2, 2, 2>(5.0e-3);
  test_rotscaletrans_control_system<3, 2, 2>(5.0e-3);
  test_rotscaletrans_control_system<2, 3, 2>();
  test_rotscaletrans_control_system<2, 2, 3>(5.0e-3);
  test_rotscaletrans_control_system<3, 3, 2>();
  test_rotscaletrans_control_system<2, 3, 3>();
  test_rotscaletrans_control_system<3, 3, 3>();
}
}  // namespace control_system
