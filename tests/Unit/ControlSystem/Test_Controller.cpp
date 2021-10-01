// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>

#include "ControlSystem/Controller.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Framework/TestHelpers.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace {
void test_controller() {
  INFO("Test controller");
  const double decrease_timescale_threshold = 1.0e-2;
  const double increase_timescale_threshold = 1.0e-4;
  const double increase_factor = 1.01;
  const double decrease_factor = 0.99;
  const double max_timescale = 10.0;
  const double min_timescale = 1.0e-3;
  const double initial_timescale = 1.0e-2;

  TimescaleTuner tst({initial_timescale}, max_timescale, min_timescale,
                     decrease_timescale_threshold, increase_timescale_threshold,
                     increase_factor, decrease_factor);

  // test following a sinusoidal target function
  double t = 0.1;
  const double dt = 1.0e-3;
  const double final_time = 5.0;
  constexpr size_t deriv_order = 2;
  const double freq = 3.0;

  // properly initialize the function of time to match our target function
  std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> f_of_t =
      std::make_unique<
          domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>>(
          t,
          std::array<DataVector, deriv_order + 1>{
              {{std::sin(freq * t)},
               {freq * std::cos(freq * t)},
               {-square(freq) * std::sin(freq * t)}}},
          t + dt);

  Controller<deriv_order> control_signal{0.5};
  const double t_offset = 0.0;

  while (t < final_time) {
    std::array<DataVector, deriv_order> target_func{
        {{std::sin(freq * t)},
         {freq * std::cos(freq * t)}}};
    const auto lambda = f_of_t->func_and_deriv(t);
    // check that the error is within the specified tolerance, which is
    // maintained by the TimescaleTuner adjusting the damping time
    CHECK(fabs(target_func[0][0] - lambda[0][0]) <
          decrease_timescale_threshold);

    // compute q = target - lambda
    // explicitly computing the derivatives of q here for testing purposes,
    // whereas, in practice, these will be computed numerically
    const auto q_and_derivs = target_func - lambda;

    // get the control signal for updating the FunctionOfTime
    const DataVector U = control_signal(t, tst.current_timescale(),
                                        q_and_derivs, t_offset, t_offset);

    t += dt;
    f_of_t->update(t, {U}, t + dt);

    // update the timescale
    tst.update_timescale({{q_and_derivs[0], q_and_derivs[1]}});
  }
}

void test_timeoffsets() {
  INFO("Test time offsets");
  const double decrease_timescale_threshold = 1.0e-2;
  const double increase_timescale_threshold = 1.0e-4;
  const double increase_factor = 1.01;
  const double decrease_factor = 0.99;
  const double max_timescale = 10.0;
  const double min_timescale = 1.0e-3;
  const double initial_timescale = 1.0e-2;

  TimescaleTuner tst({initial_timescale}, max_timescale, min_timescale,
                     decrease_timescale_threshold, increase_timescale_threshold,
                     increase_factor, decrease_factor);

  // test following a sinusoidal target function
  double t = 0.1;
  const double dt = 1.0e-3;
  const double final_time = 5.0;
  constexpr size_t deriv_order = 2;
  const double freq = 3.0;

  // some vars for a rough averaging procedure
  const double alpha = 0.1;
  std::array<DataVector, deriv_order> avg_qs{{{0.0}, {0.0}}};
  double avg_time = 0.0;

  // properly initialize the function of time to match our target function
  std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> f_of_t =
      std::make_unique<
          domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>>(
          t,
          std::array<DataVector, deriv_order + 1>{
              {{std::sin(freq * t)},
               {freq * std::cos(freq * t)},
               {-square(freq) * std::sin(freq * t)}}},
          t + dt);

  Controller<deriv_order> control_signal{0.5};

  while (t < final_time) {
    std::array<DataVector, deriv_order> target_func{
        {{std::sin(freq * t)},
         {freq * std::cos(freq * t)}}};
    const auto lambda = f_of_t->func_and_deriv(t);
    // check that the error is within the specified tolerance, which is
    // maintained by the TimescaleTuner adjusting the damping time
    CHECK(fabs(target_func[0][0] - lambda[0][0]) <
          decrease_timescale_threshold);

    // compute q = target - lambda
    // explicitly computing the derivatives of q here for testing purposes,
    // whereas, in practice, these will be computed numerically
    const auto q_and_derivs = target_func - lambda;

    // average q and its derivatives
    if (t == 0.1) {
      avg_time = t;
      avg_qs = q_and_derivs;
    } else {
      avg_time = alpha * t + (1.0 - alpha) * avg_time;
      avg_qs[0] = alpha * q_and_derivs[0] + (1.0 - alpha) * avg_qs[0];
      avg_qs[1] = alpha * q_and_derivs[1] + (1.0 - alpha) * avg_qs[1];
    }

    // get the time offset due to averaging
    const double t_offset = t - avg_time;

    // get the control signal for updating the FunctionOfTime
    const DataVector U =
        control_signal(t, tst.current_timescale(), avg_qs, t_offset, t_offset);

    t += dt;
    f_of_t->update(t, {U}, t + dt);

    // update the timescale
    tst.update_timescale({{avg_qs[0], avg_qs[1]}});
  }
}

void test_timeoffsets_noaverageq() {
  INFO("Test time offsets not averaging Q");
  const double decrease_timescale_threshold = 1.0e-2;
  const double increase_timescale_threshold = 1.0e-4;
  const double increase_factor = 1.01;
  const double decrease_factor = 0.99;
  const double max_timescale = 10.0;
  const double min_timescale = 1.0e-3;
  const double initial_timescale = 1.0e-2;

  TimescaleTuner tst({initial_timescale}, max_timescale, min_timescale,
                     decrease_timescale_threshold, increase_timescale_threshold,
                     increase_factor, decrease_factor);

  // test following a sinusoidal target function
  double t = 0.1;
  const double dt = 1.0e-3;
  const double final_time = 5.0;
  constexpr size_t deriv_order = 2;
  const double freq = 3.0;

  // some vars for a rough averaging procedure
  const double alpha = 0.1;
  std::array<DataVector, deriv_order> avg_qs{{{0.0}, {0.0}}};
  double avg_time = 0.0;

  // properly initialize the function of time to match our target function
  std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> f_of_t =
      std::make_unique<
          domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>>(
          t,
          std::array<DataVector, deriv_order + 1>{
              {{std::sin(freq * t)},
               {freq * std::cos(freq * t)},
               {-square(freq) * std::sin(freq * t)}}},
          t + dt);

  Controller<deriv_order> control_signal{0.5};

  while (t < final_time) {
    std::array<DataVector, deriv_order> target_func{
        {{std::sin(freq * t)},
         {freq * std::cos(freq * t)}}};
    const auto lambda = f_of_t->func_and_deriv(t);
    // check that the error is within the specified tolerance, which is
    // maintained by the TimescaleTuner adjusting the damping time
    CHECK(fabs(target_func[0][0] - lambda[0][0]) <
          decrease_timescale_threshold);

    // compute q = target - lambda
    // explicitly computing the derivatives of q here for testing purposes,
    // whereas, in practice, these will be computed numerically
    const auto q_and_derivs = target_func - lambda;

    // average the derivatives of q (do not average q)
    if (t == 0.1) {
      avg_time = t;
      avg_qs = q_and_derivs;
    } else {
      avg_time = alpha * t + (1.0 - alpha) * avg_time;
      avg_qs[0] = q_and_derivs[0];
      avg_qs[1] = alpha * q_and_derivs[1] + (1.0 - alpha) * avg_qs[1];
    }

    // since q is not averaged, there is no time offset
    const double q_t_offset = 0.0;
    // get the derivative time offset due to averaging
    const double t_offset = t - avg_time;

    // get the control signal for updating the FunctionOfTime
    const DataVector U = control_signal(t, tst.current_timescale(), avg_qs,
                                        q_t_offset, t_offset);

    t += dt;
    f_of_t->update(t, {U}, t + dt);

    // update the timescale
    tst.update_timescale({{avg_qs[0], avg_qs[1]}});
  }
}

void test_equality_and_serialization() {
  INFO("Test equality and serialization");
  Controller<2> controller1{0.5};
  Controller<2> controller2{1.0};
  Controller<2> controller3{0.5};
  controller3.assign_time_between_updates(2.0);

  CHECK(controller1 != controller2);
  CHECK_FALSE(controller1 == controller2);
  CHECK(controller1 != controller3);

  controller1.assign_time_between_updates(2.0);

  CHECK(controller1 == controller3);
  CHECK(controller1.get_update_fraction() == 0.5);

  Controller<2> controller1_serialized = serialize_and_deserialize(controller1);

  CHECK(controller1 == controller1_serialized);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.Controller", "[ControlSystem][Unit]") {
  test_controller();
  test_timeoffsets();
  test_timeoffsets_noaverageq();
  test_equality_and_serialization();
}
