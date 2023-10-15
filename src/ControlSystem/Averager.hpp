// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <deque>
#include <optional>
#include <pup.h>
#include <pup_stl.h>

#include "DataStructures/DataVector.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"

/// \ingroup ControlSystemGroup
/// A weighted exponential averager of \f$Q\f$ and its derivatives
/// implementing Appendix A in \cite Hemberger2012jz
///
/// The purpose of the averager is to provide \f$Q\f$ and 'smoothed' numerical
/// derivatives of \f$Q\f$ up to the `DerivOrder`'th derivative. The 0th
/// derivative of \f$Q\f$ is not typically averaged, since no differencing is
/// necessary and no noise is introduced. The `average_0th_deriv_of_q` option
/// allows for specifying that the 0th derivative should be averaged in addition
/// to the derivatives. This may be desirable for control systems where
/// \f$Q\f$ contains some 'noise', e.g. size control typically uses an averaged
/// \f$Q\f$ since the raw computed \f$Q\f$ is a function of the minimum on a
/// surface and may jump around discontinuously.
///
/// The averager is designed to support DerivOrders 1, 2, and 3. If an
/// additional DerivOrder is needed, finite differencing needs to be implemented
/// to specifically handle that order (as it seems like overkill to generalize
/// the differencing stencil at this time).
template <size_t DerivOrder>
class Averager {
  static_assert(DerivOrder != 0, "Cannot average just a function value.");
  static_assert(
      DerivOrder <= 3,
      "Averager only instantiated for deriv order less than or equal to 3. If "
      "you want a higher deriv order, you'll need to add it yourself.");

 public:
  struct AverageTimescaleFraction {
    using type = double;
    inline const static std::string help {
        "Time scale of exponential averaging"};
  };

  struct Average0thDeriv {
    using type = bool;
    inline const static std::string help {
        "Whether to average the 0th derivative"};
  };

  using options = tmpl::list<AverageTimescaleFraction, Average0thDeriv>;
  inline const static std::string help{
      "Averager: Performs exponential averaging of the control signal at "
      "multiple times in order to provide smoother derivatives of the control "
      "signal."};

  /// `avg_timescale_frac` determines the exponential averaging timescale
  /// through \f$\tau_\mathrm{avg} = \f$`avg_timescale_frac`\f$\times \tau\f$,
  /// where \f$\tau\f$ is the damping time.
  /// `avg_timescale_frac` must be positive.
  /// `average_0th_deriv_of_q` determines whether the call operator returns an
  /// averaged or unaveraged quantity for the 0th derivative of \f$Q\f$.
  /// `true` returns an averaged 0th derivative of \f$Q\f$.
  /// `false` returns the raw 0th derivative of \f$Q\f$.
  /// The derivatives are always averaged (to reduce noise due to numerical
  /// differentiation), so the `average_0th_deriv_of_q` option only specifies
  /// whether to return an averaged value for the 0th derivative piece of the
  /// function.
  Averager(double avg_timescale_frac, bool average_0th_deriv_of_q);

  Averager() = default;
  // Explicitly defined move constructor due to the fact that the std::deque
  // move constructor is not marked
  Averager(Averager&& rhs);
  Averager& operator=(Averager&& rhs);
  Averager(const Averager&) = default;
  Averager& operator=(const Averager&) = default;
  ~Averager() = default;

  template <size_t DDerivOrder>                          // NOLINT
  friend bool operator==(const Averager<DDerivOrder>&,   // NOLINT
                         const Averager<DDerivOrder>&);  // NOLINT

  /// Returns \f$Q\f$ and its derivatives at \f$t=\f$`time`, provided there is
  /// sufficient data. The averager is limited by the need for at least
  /// (`DerivOrder` + 1) data points in order to provide the `DerivOrder`'th
  /// derivative. If sufficient data is available, it returns \f$Q\f$ and
  /// averaged derivatives of \f$Q\f$ up to the `DerivOrder`'th derivative, at
  /// \f$t=\f$`time`. If `using_average_0th_deriv_of_q()` is `true`, then the
  /// returned 0th derivative of \f$Q\f$ is also averaged. In the case that
  /// there is insufficient data, the operator returns an invalid
  /// `std::optional`.
  const std::optional<std::array<DataVector, DerivOrder + 1>>& operator()(
      double time) const;
  /// A function that allows for resetting the averager.
  void clear();
  /// The function responsible for updating the averaged values
  /// at \f$t=\f$`time`. Requires `raw_q` (the raw components of \f$Q(t)\f$)
  /// and `timescales` (the associated damping times for each component).
  void update(double time, const DataVector& raw_q,
              const DataVector& timescales);
  /// Returns the latest time at which the averager has sufficient data to
  /// return \f$Q\f$ and its derivatives.
  double last_time_updated() const;
  /// Returns the exponentially averaged time at \f$t=\f$`time`. The time is
  /// averaged along side \f$Q\f$ to determine the effective time at which
  /// the average value is computed. The effective time is retarded, due to the
  /// weighting of past times.
  double average_time(double time) const;
  /// Returns a bool corresponding to whether `average_0th_deriv_of_q`
  /// is `true`/`false`.
  bool using_average_0th_deriv_of_q() const { return average_0th_deriv_of_q_; };
  /// Returns the averaging timescale fraction
  double avg_timescale_frac() const { return avg_tscale_frac_; }

  /// Assign the minimum of the measurement timescales related to a specific
  /// control system to time_between_measurements_
  void assign_time_between_measurements(
      const double current_time_between_measurements) {
    time_between_measurements_ = current_time_between_measurements;
  }

  /// Returns `true` if the averager is ready to receive a measurement
  bool is_ready(const double time) const {
    if (UNLIKELY(times_.empty())) {
      return true;
    }
    return time >= times_[0] + time_between_measurements_;
  }

  void pup(PUP::er& p);

 private:
  /// Returns the function and numerical derivatives up to the
  /// `DerivOrder`'th derivative using backward finite differencing that
  /// supports non-uniform time spacing. Since we use a stencil size of
  /// (`DerivOrder` + 1), the order of the finite difference varies for
  /// different derivatives. Assuming \f$M\lte\f$`DerivOrder`, the \f$M\f$'th
  /// derivative is of order (`DerivOrder` + 1 - \f$M\f$).
  std::array<DataVector, DerivOrder + 1> get_derivs() const;

  double avg_tscale_frac_;
  bool average_0th_deriv_of_q_;
  std::optional<std::array<DataVector, DerivOrder + 1>> averaged_values_{};
  std::optional<std::array<DataVector, DerivOrder + 1>> nullopt_ = std::nullopt;
  std::deque<double> times_;
  std::deque<DataVector> raw_qs_;
  double weight_k_ = 0.0;
  double tau_k_ = 0.0;
  // If this time_between_measurements_ isn't set, the default should just be
  // that no measurements happen (i.e. infinity)
  double time_between_measurements_{std::numeric_limits<double>::infinity()};
};

template <size_t DDerivOrder>
bool operator!=(const Averager<DDerivOrder>&, const Averager<DDerivOrder>&);
