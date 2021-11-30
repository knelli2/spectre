// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "ApparentHorizons/Strahlkorper.hpp"
#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Shape.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Systems/Shape.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Shape.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/ControlSystem/SystemHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Frame {
struct Distorted;
}  // namespace Frame

namespace control_system {
namespace {
template <size_t DerivOrder>
struct MockMetavars {
  static constexpr size_t volume_dim = 3;

  enum class Phase { Initialization, Testing, Exit };

  using metavars = MockMetavars<DerivOrder>;

  using shape_system =
      control_system::Systems::Shape<control_system::ah::HorizonLabel::AhA,
                                     DerivOrder>;

  using element_component = TestHelpers::MockElementComponent<metavars>;
  using shape_component =
      TestHelpers::MockControlComponent<metavars, shape_system>;

  using component_list = tmpl::list<element_component, shape_component>;
};

using FoTMap = std::unordered_map<
    std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;
// We use the grid frame here because that is where the apparent horizons will
// be found, however, for this test it doesn't really matter because we
// randomize the initial data
using Strahlkorper = Strahlkorper<Frame::Distorted>;

template <size_t DerivOrder>
void test_shape_control(
    const double initial_time, const double final_time, Strahlkorper& fake_ah,
    const domain::FunctionsOfTime::PiecewisePolynomial<DerivOrder>&
        ah_coefs_function_of_time,
    gsl::not_null<std::mt19937*> generator, const double func_eps,
    const double deriv_eps) {
  using metavars = MockMetavars<DerivOrder>;
  using system = typename metavars::shape_system;
  using shape_component = typename metavars::shape_component;
  using element_component = typename metavars::element_component;
  using shape_init_simple_tags = TestHelpers::init_simple_tags<system>;
  auto& mutable_ah_coefs = fake_ah.coefficients();
  const auto& fake_ah_coefs = fake_ah.coefficients();

  // Setup initial shape map coefficients. In the map the coefficients are
  // stored as the negative of the actual spherical harmonic coefficients
  // because that's just how the map is defined. But since these are random
  // numbers it doesn't matter for initial data
  auto initial_shape_func = make_array<DerivOrder + 1, DataVector>(
      DataVector{fake_ah_coefs.size(), 0.0});
  SpherepackIterator iter{fake_ah.l_max(), fake_ah.m_max()};
  std::uniform_real_distribution<double> coef_dist{-1.0, 1.0};
  for (size_t i = 0; i < initial_shape_func.size(); i++) {
    for (iter.reset(); iter; ++iter) {
      // Enforce l=0,l=1 components to be 0 always
      if (iter.l() == 0 or iter.l() == 1) {
        continue;
      }
      gsl::at(initial_shape_func, i)[iter()] =
          make_with_random_values<double>(generator, coef_dist, 1);
    }
  }

  // Setup initial size map parameters
  auto initial_size_func =
      make_array<DerivOrder + 1, DataVector>(DataVector{1, 0.0});
  const double ah_radius = fake_ah.average_radius();
  // Excision sphere radius needs to be inside the AH
  const double excision_radius = 0.75 * ah_radius;
  // We only test for a constant size function of time. A test with a changing
  // size can be added in later
  initial_size_func[0][0] = ah_radius;

  // Setup control system stuff
  const std::string shape_name =
      Systems::Shape<ah::HorizonLabel::AhA, DerivOrder>::name();
  const std::string size_name =
      ControlErrors::detail::size_name<ah::HorizonLabel::AhA>();
  const std::string excision_sphere_A_name =
      ControlErrors::detail::excision_sphere_name<ah::HorizonLabel::AhA>();
  Averager<DerivOrder - 1> shape_averager{0.25, false};
  Controller<DerivOrder> shape_controller{0.3};
  const double update_timescale = 0.5;
  TimescaleTuner shape_tuner{
      std::vector<double>(fake_ah_coefs.size(), update_timescale),
      10.0,
      0.1,
      4.0e-3,
      1.0e-3,
      1.01,
      0.98};
  ControlErrors::Shape<ah::HorizonLabel::AhA> shape_control_error{};
  const double shape_min_timescale = min(shape_tuner.current_timescale());
  shape_controller.assign_time_between_updates(shape_min_timescale);
  shape_controller.set_initial_time(initial_time);
  const std::array<DataVector, 1> shape_measurement_timescale{
      {control_system::Tags::detail::calculate_measurement_timescales(
          shape_controller, shape_tuner)}};
  shape_averager.assign_time_between_measurements(
      min(shape_measurement_timescale[0]));
  LinkedMessageQueue<double,
                     tmpl::list<QueueTags::Strahlkorper<Frame::Distorted>>>
      empty_queue{};
  auto init_shape_tuple =
      tuples::tagged_tuple_from_typelist<shape_init_simple_tags>{
          shape_averager,   shape_tuner,         shape_name,
          shape_controller, shape_control_error, empty_queue};
  const double initial_shape_expiration_time =
      shape_controller.get_update_fraction() * shape_min_timescale;

  // Since the map for A/B are independent of each other, we only need to test
  // one of them
  FoTMap initial_functions_of_time{};
  FoTMap initial_measurement_timescales{};
  initial_functions_of_time[shape_name] = std::make_unique<
      domain::FunctionsOfTime::PiecewisePolynomial<DerivOrder>>(
      initial_time, initial_shape_func, initial_shape_expiration_time);
  initial_functions_of_time[size_name] = std::make_unique<
      domain::FunctionsOfTime::PiecewisePolynomial<DerivOrder>>(
      initial_time, initial_size_func, std::numeric_limits<double>::infinity());
  initial_measurement_timescales[shape_name] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          initial_time, shape_measurement_timescale,
          initial_shape_expiration_time);

  // Fake domain
  Domain<3> fake_domain{
      {},
      {},
      {{excision_sphere_A_name,
        ExcisionSphere<3>{excision_radius, fake_ah.center()}}}};

  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavars>;
  MockRuntimeSystem runner{{std::move(fake_domain)},
                           {std::move(initial_functions_of_time),
                            std::move(initial_measurement_timescales)}};
  ActionTesting::emplace_singleton_component_and_initialize<shape_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, init_shape_tuple);
  ActionTesting::emplace_array_component<element_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0);

  ActionTesting::set_phase(make_not_null(&runner), metavars::Phase::Testing);

  auto& cache = ActionTesting::cache<element_component>(runner, 0);
  const auto& functions_of_time =
      Parallel::get<domain::Tags::FunctionsOfTime>(cache);
  const auto& measurement_timescales =
      Parallel::get<control_system::Tags::MeasurementTimescales>(cache);

  double time = initial_time;
  std::optional<double> prev_time{};
  LinkedMessageId<double> measurement_id{};
  while (time < final_time) {
    // Setup the measurement id. This would normally be created in the control
    // system event.
    measurement_id = LinkedMessageId<double>{time, prev_time};

    // Update expected strahlkorper coefficients. This serves as our measurement
    // in the distorted frame
    mutable_ah_coefs = ah_coefs_function_of_time.func(time)[0];

    // Only need one measurement technically but in practice we will receive
    // two (one for each horizon). The measurement for the wrong horizon
    // shouldn't do anything so we just give it the same strahlkorper
    system::process_measurement::apply(
        ah::BothHorizons::FindHorizon<ah::HorizonLabel::AhA>{}, fake_ah, cache,
        measurement_id);
    system::process_measurement::apply(
        ah::BothHorizons::FindHorizon<ah::HorizonLabel::AhB>{}, fake_ah, cache,
        measurement_id);
    // Only one action should be queued (the one for A)
    CHECK(ActionTesting::number_of_queued_simple_actions<shape_component>(
              runner, 0) == 1);
    ActionTesting::invoke_queued_simple_action<shape_component>(
        make_not_null(&runner), 0);

    // Our dt is set by the smallest measurement timescale. The control system
    // updates these timescales when it updates the functions of time
    prev_time = time;
    time += min(measurement_timescales.at(shape_name)->func(time)[0]);
  }

  const auto lambda_00_coef =
      functions_of_time.at(size_name)->func(final_time)[0][0];
  const double Y00_coef = ControlErrors::detail::y00_coef();

  auto ah_coefs_and_derivs =
      ah_coefs_function_of_time.func_and_2_derivs(final_time);

  // Our expected coefs are just the (minus) coefs of the AH (except for l=0,l=1
  // which should be zero) scaled by the relative size factor defined in the
  // control error
  auto expected_shape_coefs =
      -1.0 * (excision_radius / Y00_coef - lambda_00_coef) /
      ah_coefs_and_derivs[0][iter.set(0, 0)()] * ah_coefs_and_derivs;
  // Manually set 0,0 component to 0
  expected_shape_coefs[0][iter.set(0, 0)()] = 0.0;

  const auto& shape_func =
      dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<DerivOrder>&>(
          *functions_of_time.at(shape_name));

  const auto shape_coefs = shape_func.func_and_2_derivs(final_time);

  Approx custom_approx = Approx::custom().epsilon(1.0).scale(1.0);
  for (size_t i = 0; i < shape_coefs.size(); i++) {
    if (i == 0) {
      custom_approx = Approx::custom().epsilon(func_eps).scale(1.0);
    } else {
      custom_approx = Approx::custom().epsilon(deriv_eps).scale(1.0);
    }
    INFO("i = " + get_output(i));
    CHECK_ITERABLE_CUSTOM_APPROX(gsl::at(shape_coefs, i),
                                 gsl::at(expected_shape_coefs, i),
                                 custom_approx);
  }
}

template <size_t DerivOrder, typename Generator>
void test_suite(const double initial_time, const double final_time,
                const double radius, Strahlkorper& fake_ah,
                gsl::not_null<Generator*> generator, const double func_eps,
                const double deriv_eps) {
  SpherepackIterator iter{fake_ah.l_max(), fake_ah.m_max()};
  std::uniform_real_distribution<double> coef_dist{-1.0, 1.0};
  const DataVector zero_dv{fake_ah.coefficients().size(), 0.0};
  const auto zero_array = make_array<DerivOrder + 1, DataVector>(zero_dv);

  // Allocate now and reuse as we go
  auto initial_ah_coefs = zero_array;
  domain::FunctionsOfTime::PiecewisePolynomial<DerivOrder> expected_ah_coefs{};
  auto& fake_ah_coefs = fake_ah.coefficients();

  const std::string deriv_order_string =
      "DerivOrder=" + get_output(DerivOrder) + ": ";
  {
    INFO(deriv_order_string + "Stationary spherical AH");
    // Set all coefs back to zero first
    initial_ah_coefs = zero_array;
    initial_ah_coefs[0][iter.set(0, 0)()] = radius;
    fake_ah_coefs = initial_ah_coefs[0];
    expected_ah_coefs = {initial_time, initial_ah_coefs,
                         std::numeric_limits<double>::infinity()};
    test_shape_control<DerivOrder>(initial_time, final_time, fake_ah,
                                   expected_ah_coefs, generator, deriv_eps,
                                   deriv_eps);
  }
  {
    INFO(deriv_order_string + "Stationary non-spherical AH");
    // Set all coefs back to zero first
    initial_ah_coefs = zero_array;
    for (iter.reset(); iter; ++iter) {
      // Enforce l=0,l=1 components to be 0 always
      if (iter.l() == 0 or iter.l() == 1) {
        continue;
      }
      initial_ah_coefs[0][iter()] =
          make_with_random_values<double>(generator, coef_dist, 1);
    }
    // Ensure that radius of AH is positive and constant
    initial_ah_coefs[0][iter.set(0, 0)()] = radius;
    initial_ah_coefs[1] = zero_dv;
    initial_ah_coefs[2] = zero_dv;
    fake_ah_coefs = initial_ah_coefs[0];
    expected_ah_coefs = {initial_time, initial_ah_coefs,
                         std::numeric_limits<double>::infinity()};
    test_shape_control<DerivOrder>(initial_time, final_time, fake_ah,
                                   expected_ah_coefs, generator, deriv_eps,
                                   deriv_eps);
    (void)func_eps;
  }
  {
    INFO(deriv_order_string + "AH coefficients linearly increasing/decreasing");
    // Set all coefs back to zero first
    initial_ah_coefs = zero_array;
    for (iter.reset(); iter; ++iter) {
      // Enforce l=0,l=1 components to be 0 always
      if (iter.l() == 0 or iter.l() == 1) {
        continue;
      }
      initial_ah_coefs[0][iter()] =
          make_with_random_values<double>(generator, coef_dist, 1);
      initial_ah_coefs[1][iter()] =
          make_with_random_values<double>(generator, coef_dist, 1);
    }
    // Ensure that radius of AH is positive and constant
    initial_ah_coefs[0][iter.set(0, 0)()] = radius;
    initial_ah_coefs[1][iter.set(0, 0)()] = 0.0;
    initial_ah_coefs[2] = zero_dv;
    fake_ah_coefs = initial_ah_coefs[0];
    expected_ah_coefs = {initial_time, initial_ah_coefs,
                         std::numeric_limits<double>::infinity()};
    test_shape_control<DerivOrder>(initial_time, final_time, fake_ah,
                                   expected_ah_coefs, generator, deriv_eps,
                                   deriv_eps);
  }
  {
    INFO(deriv_order_string +
         "AH coefficients quadratically increasing/decreasing");
    // Set all coefs back to zero first
    initial_ah_coefs = zero_array;
    for (iter.reset(); iter; ++iter) {
      // Enforce l=0,l=1 components to be 0 always
      if (iter.l() == 0 or iter.l() == 1) {
        continue;
      }
      for (size_t i = 0; i < initial_ah_coefs.size(); i++) {
        gsl::at(initial_ah_coefs, i)[iter()] =
            make_with_random_values<double>(generator, coef_dist, 1);
        if (i == initial_ah_coefs.size() - 1) {
          gsl::at(initial_ah_coefs, i)[iter()] /= cube(double(i + 5));
        }
      }
    }
    // Ensure that radius of AH is positive and constant
    initial_ah_coefs[0][iter.set(0, 0)()] = radius;
    initial_ah_coefs[1][iter.set(0, 0)()] = 0.0;
    initial_ah_coefs[2][iter.set(0, 0)()] = 0.0;
    fake_ah_coefs = initial_ah_coefs[0];
    expected_ah_coefs = {initial_time, initial_ah_coefs,
                         std::numeric_limits<double>::infinity()};
    test_shape_control<DerivOrder>(initial_time, final_time, fake_ah,
                                   expected_ah_coefs, generator, func_eps,
                                   deriv_eps);
  }
  {
    INFO(deriv_order_string + "AH coefficients oscillating sinusoidally");
    // Set all coefs back to zero first
    initial_ah_coefs = zero_array;
    // All coefficients (and derivatives) are of the form
    // coef = sin(freq * time + offset)
    // dtcoef = freq * cos(freq * time + offset)
    // d2tcoef = -freq^2 * sin(freq * time + offset)
    // d3tcoef = -freq^3 * cos(freq * time + offset)

    // Only allow at most one full oscillation by the time we reach the end.
    // That way the coefficients aren't oscillating too fast
    std::uniform_real_distribution<double> freq_dist{0.0,
                                                     2.0 * M_PI / final_time};
    double offset = 0.0;
    double freq = 0.0;
    for (iter.reset(); iter; ++iter) {
      // Enforce l=0,l=1 components to be 0 always
      if (iter.l() == 0 or iter.l() == 1) {
        continue;
      }
      offset = make_with_random_values<double>(generator, coef_dist, 1);
      freq = make_with_random_values<double>(generator, freq_dist, 1);
      initial_ah_coefs[0][iter()] = sin(freq * initial_time + offset);
      initial_ah_coefs[1][iter()] = freq * cos(freq * initial_time + offset);
      initial_ah_coefs[2][iter()] =
          -square(freq) * sin(freq * initial_time + offset);
      if (DerivOrder > 2) {
        initial_ah_coefs[3][iter()] =
            -cube(freq) * cos(freq * initial_time + offset);
      }
    }
    // Ensure that radius of AH is positive and constant
    initial_ah_coefs[0][iter.set(0, 0)()] = radius;
    initial_ah_coefs[1][iter.set(0, 0)()] = 0.0;
    initial_ah_coefs[2][iter.set(0, 0)()] = 0.0;
    fake_ah_coefs = initial_ah_coefs[0];

    // Initialize the expected function of time
    const double dt = 0.01;
    double time = initial_time + dt;
    expected_ah_coefs = {initial_time, initial_ah_coefs, time};

    // Update it's derivative often so the function is smooth
    DataVector updated_deriv = zero_dv;
    for (iter.reset(); iter; ++iter) {
      if (iter.l() == 0 or iter.l() == 1) {
        continue;
      }
      if (DerivOrder == 2) {
        updated_deriv[iter()] = -square(freq) * sin(freq * time + offset);
      } else {
        updated_deriv[iter()] = -cube(freq) * cos(freq * time + offset);
      }
    }
    while (time < final_time) {
      expected_ah_coefs.update(time, updated_deriv, time + dt);
      time += dt;
    }

    test_shape_control<DerivOrder>(initial_time, final_time, fake_ah,
                                   expected_ah_coefs, generator, func_eps,
                                   deriv_eps);
  }
}

SPECTRE_TEST_CASE("Unit.ControlSystem.Systems.Shape", "[ControlSystem][Unit]") {
  MAKE_GENERATOR(generator);
  domain::FunctionsOfTime::register_derived_with_charm();

  const std::array<double, 3> origin{{0.0, 0.0, 0.0}};
  const double radius = 1.0;
  const double initial_time = 0.0;
  const double final_time = 150.0;
  // For some of the faster changing AHs, (quadratic and sinusoid) the control
  // system isn't able to perfectly match the expected map parameters, but it is
  // able to get the control error to be a small constant offset (rather than
  // 0). This larger epsilon accounts for the difference in the map parameters
  // (not their derivatives though which it is able to more accurately match)
  // based off the constant offset in the control error.
  const double custom_approx_func_eps = 5.0e-3;
  // This epsilon is mostly necessary for deriv order 3. For most of deriv order
  // 2, we can use an epsilon of 5.0e-10, but occasionally there are cases where
  // the error for deriv order 2 is larger than 5.0e-10 so to be safe we just
  // use 5.0e-5 all the time.
  const double custom_approx_deriv_eps = 5.0e-5;

  std::vector<size_t> ells_to_test{2, 3, 5, 10, 12, 13};
  // std::vector<size_t> ells_to_test{2};
  for (const auto& ell : ells_to_test) {
    Strahlkorper fake_ah{ell, ell, radius, origin};

    test_suite<2>(initial_time, final_time, radius, fake_ah,
                  make_not_null(&generator), custom_approx_func_eps,
                  custom_approx_deriv_eps);
    test_suite<3>(initial_time, final_time, radius, fake_ah,
                  make_not_null(&generator), custom_approx_func_eps,
                  custom_approx_deriv_eps);
  }
}
}  // namespace
}  // namespace control_system
