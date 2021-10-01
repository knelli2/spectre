// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <ostream>

#include "DataStructures/DataVector.hpp"
#include "Options/Options.hpp"
#include "Parallel/Printf.hpp"

/// \cond
namespace PUP {
class er;
}
/// \endcond

/// \ingroup ControlSystemGroup
/// A PND (proportional to Q and N derivatives of Q) controller that computes
/// the control signal:
/// \f[ U(t) = \sum_{k=0}^{N} a_{k} \frac{d^kQ}{dt^k} \f]
/// where N is specified by the template parameter `DerivOrder`.
///
/// If an averager is used for `q_and_derivs` (as we typically do), there is an
/// induced time offset, \f$\Delta t\f$, due to the time-weighted averaging.
/// Therefore, the `q_and_derivs` that we have in hand are at some time
/// \f$t_{0}\f$. However, we desire `q_and_derivs` at the current time
/// \f$t = t_{0} + \Delta t\f$ to determine the appropriate control
/// signal. We accomplish this by Taylor expanding
/// \f$Q(t_{0} + \Delta t)\f$. The averager allows for averaging of
/// \f$Q\f$ and its derivatives OR to not average \f$Q\f$ while still averaging
/// the derivatives (the derivatives are always averaged in order to reduce
/// noise due to numerical differentiation). When they are both averaged, the
/// time offset will be identical for \f$Q\f$ and the derivatives,
/// i.e. `q_time_offset` = `deriv_time_offset`. If an unaveraged \f$Q\f$ is
/// used, then the time offset associated with \f$Q\f$ is zero,
/// i.e. `q_time_offset`=0. and the derivative time offset, `deriv_time_offset`,
/// remains non-zero.
template <size_t DerivOrder>
class Controller {
 public:
  struct UpdateFraction {
    using type = double;
    static constexpr Options::String help = {
        "Fraction of damping timescale used to update functions of time."};

    static type suggested_value() noexcept { return 0.3; }
  };

  using options = tmpl::list<UpdateFraction>;
  static constexpr Options::String help{
      "Controller. Computes control signal used to reset highest derivative of "
      "a function of time."};

  Controller(const double update_fraction) noexcept
      : update_fraction_(update_fraction) {}

  Controller() = default;
  Controller(Controller&&) noexcept = default;
  Controller& operator=(Controller&&) noexcept = default;
  Controller(const Controller&) = default;
  Controller& operator=(const Controller&) = default;
  ~Controller() = default;

  DataVector operator()(
      const double time, const DataVector& timescales,
      const std::array<DataVector, DerivOrder + 1>& q_and_derivs,
      double q_time_offset, double deriv_time_offset) const noexcept;

  bool is_triggered(const double time) const noexcept {
    return time >= last_trigger_ + time_between_triggers_;
  }

  void assign_time_between_triggers(
      const double current_timescale) const noexcept {
    time_between_triggers_ = update_fraction_ * current_timescale;
  }

  double next_expiration_time(
      const double current_expiration_time) const noexcept {
    return current_expiration_time + time_between_triggers_;
  }

  std::string get_output() const noexcept {
    std::ostringstream os;
    os << "last_trigger_ = " << last_trigger_ << "\n";
    os << "update_fraction_ = " << update_fraction_ << "\n";
    os << "time_between_triggers_ = " << time_between_triggers_ << "\n";
    return os.str();
  }

  void pup(PUP::er& p) {
    p | time_between_triggers_;
    p | last_trigger_;
  }

 private:
  // double
  // time_between_triggers_{std::numeric_limits<double>::signaling_NaN()};
  mutable double time_between_triggers_;
  mutable double last_trigger_{0.0};
  // double update_fraction_{std::numeric_limits<double>::signaling_NaN()};
  double update_fraction_{0.3};
};
