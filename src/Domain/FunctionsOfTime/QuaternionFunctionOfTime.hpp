// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <boost/math/quaternion.hpp>
#include <cmath>
#include <limits>
#include <pup.h>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTimeHelpers.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Utilities/Gsl.hpp"

namespace domain::FunctionsOfTime {
/// \ingroup ComputationalDomainGroup
/// \brief A FunctionOfTime that stores quaternions for the rotation map
///
/// \details This FunctionOfTime stores quaternions that will be used in the
/// time-dependent rotation map as well as the orbital angular velocity that
/// will be controlled by the rotation control sytem. To get the quaternion, an
/// ODE is solved of the form \f$ \dot{q} = \frac{1}{2} q \times \omega \f$
/// where \f$ \omega \f$ is the orbital angular velocity which is stored
/// internally as the derivative of an angle `PiecewisePolynomial`, and
/// \f$ \times \f$ here is quaternion multiplication.
///
/// Different from a `PiecewisePolynomial`, only the quaternion
/// itself is stored, not any of the derivatives because the derivatives must be
/// calculated from the solved ODE at every function call. Because
/// derivatives of the quaternion are not stored, the template parameter
/// `MaxDeriv` refers to both the max derivative of the stored angle
/// PiecewisePolynomial and the max derivative returned by the
/// QuaternionFunctionOfTime. The `update` function is then just a wrapper
/// around the internal `PiecewisePolynomial::update` function with the addition
/// that it then updates the stored quaternions as well.
///
/// The angle PiecewisePolynomial is accessible through the `angle_func`,
/// `angle_func_and_deriv`, and `angle_func_and_2_derivs` functions which
/// correspond to the function calls of a normal PiecewisePolynomial except
/// without the `angle_` prefix.
///
/// It is encouraged to use `quat_func` and `angle_func` when you want the
/// specific values of the functions to avoid ambiguity in what you are
/// calling. However, the original three `func` functions inherited from the
/// FunctionOfTime base class are necessary because the maps use the generic
/// `func` functions, thus they return the quaternion and its derivatives (which
/// are needed for the map). This is all to keep the symmetry of naming
/// `angle_func` and `quat_func` so that function calls won't be ambiguous.
template <size_t MaxDeriv>
class QuaternionFunctionOfTime : public FunctionOfTime {
 public:
  QuaternionFunctionOfTime() = default;
  QuaternionFunctionOfTime(
      double t, std::array<DataVector, 1> initial_quat_func,
      std::array<DataVector, MaxDeriv + 1> initial_angle_func,
      double expiration_time);

  ~QuaternionFunctionOfTime() override = default;
  QuaternionFunctionOfTime(QuaternionFunctionOfTime&&) = default;
  QuaternionFunctionOfTime& operator=(QuaternionFunctionOfTime&&) = default;
  QuaternionFunctionOfTime(const QuaternionFunctionOfTime&) = default;
  QuaternionFunctionOfTime& operator=(const QuaternionFunctionOfTime&) =
      default;

  // LCOV_EXCL_START
  explicit QuaternionFunctionOfTime(CkMigrateMessage* /*unused*/) {}
  // LCOV_EXCL_STOP

  auto get_clone() const -> std::unique_ptr<FunctionOfTime> override;

  // clang-tidy: google-runtime-references
  // clang-tidy: cppcoreguidelines-owning-memory,-warnings-as-errors
  WRAPPED_PUPable_decl_template(QuaternionFunctionOfTime<MaxDeriv>);  // NOLINT

  void reset_expiration_time(const double next_expiration_time) override {
    angle_f_of_t_.reset_expiration_time(next_expiration_time);
  }

  /// Returns domain of validity for the function of time
  std::array<double, 2> time_bounds() const override {
    return angle_f_of_t_.time_bounds();
  }

  /// Updates the `MaxDeriv`th derivative of the angle piecewisepolynomial at
  /// the given time, then updates the stored quaternions.
  ///
  /// `updated_max_deriv` is a datavector of the `MaxDeriv`s for each component.
  /// `next_expiration_time` is the next expiration time.
  void update(double time_of_update, DataVector updated_max_deriv,
              double next_expiration_time) override;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

  /// Returns the quaternion at an arbitrary time `t`.
  std::array<DataVector, 1> func(const double t) const override {
    return quat_func(t);
  }

  /// Returns the quaternion and its first derivative at an arbitrary time `t`.
  std::array<DataVector, 2> func_and_deriv(const double t) const override {
    return quat_func_and_deriv(t);
  }

  /// Returns the quaternion and the first two derivatives at an arbitrary
  /// time `t`.
  std::array<DataVector, 3> func_and_2_derivs(const double t) const override {
    return quat_func_and_2_derivs(t);
  }

  /// Returns the quaternion at an arbitrary time `t`.
  std::array<DataVector, 1> quat_func(double t) const;

  /// Returns the quaternion and its first derivative at an arbitrary time `t`.
  std::array<DataVector, 2> quat_func_and_deriv(double t) const;

  /// Returns the quaternion and the first two derivatives at an arbitrary
  /// time `t`.
  std::array<DataVector, 3> quat_func_and_2_derivs(double t) const;

  /// Returns stored angle at an arbitrary time `t`.
  std::array<DataVector, 1> angle_func(const double t) const {
    return angle_f_of_t_.func(t);
  }

  /// Returns stored angle and its first derivative (omega) at an arbitrary time
  /// `t`.
  std::array<DataVector, 2> angle_func_and_deriv(const double t) const {
    return angle_f_of_t_.func_and_deriv(t);
  }

  /// Returns stored angle and the first two derivatives at an arbitrary
  /// time `t`.
  std::array<DataVector, 3> angle_func_and_2_derivs(const double t) const {
    return angle_f_of_t_.func_and_2_derivs(t);
  }

  std::array<DataVector, MaxDeriv + 1> angle_func_and_all_derivs(
      const double t) const {
    return angle_f_of_t_.func_and_all_derivs(t);
  }

 private:
  template <size_t LocalMaxDeriv>
  friend bool operator==(  // NOLINT(readability-redundant-declaration)
      const QuaternionFunctionOfTime<LocalMaxDeriv>& lhs,
      const QuaternionFunctionOfTime<LocalMaxDeriv>& rhs);

  template <size_t LocalMaxDeriv>
  friend std::ostream& operator<<(  // NOLINT(readability-redundant-declaration)
      std::ostream& os,
      const QuaternionFunctionOfTime<LocalMaxDeriv>& quaternion_f_of_t);

  std::vector<FunctionOfTimeHelpers::StoredInfo<1, false>>
      stored_quaternions_and_times_;

  domain::FunctionsOfTime::PiecewisePolynomial<MaxDeriv> angle_f_of_t_;

  /// Integrates the ODE \f$ \dot{q} = \frac{1}{2} q \times \omega \f$ from time
  /// `t0` to time `t`. On input, `quaternion_to_integrate` is the initial
  /// quaternion at time `t0` and on output, it stores the result at time `t`
  void solve_quaternion_ode(
      gsl::not_null<boost::math::quaternion<double>*> quaternion_to_integrate,
      double t0, double t) const;

  /// Updates the `std::vector<StoredInfo>` to have the same number of stored
  /// quaternions as the `angle_f_of_t_ptr` has stored angles. This is necessary
  /// to ensure we can solve the ODE at any time `t`
  void update_stored_info();

  /// Does common operations to all the `func` functions such as updating stored
  /// info, solving the ODE, and returning the normalized quaternion as a boost
  /// quaternion for easy calculations
  boost::math::quaternion<double> setup_func(double t) const;
};

template <size_t MaxDeriv>
bool operator!=(const QuaternionFunctionOfTime<MaxDeriv>& lhs,
                const QuaternionFunctionOfTime<MaxDeriv>& rhs);

/// \cond
template <size_t MaxDeriv>
PUP::able::PUP_ID QuaternionFunctionOfTime<MaxDeriv>::my_PUP_ID = 0;  // NOLINT
/// \endcond
}  // namespace domain::FunctionsOfTime
