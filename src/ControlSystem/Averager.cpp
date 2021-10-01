// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/Averager.hpp"

#include <ostream>

#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/StdArrayHelpers.hpp"

#include "Parallel/Printf.hpp"

template <size_t DerivOrder>
Averager<DerivOrder>::Averager(const double avg_timescale_frac,
                               const bool average_0th_deriv_of_q) noexcept
    : avg_tscale_frac_(avg_timescale_frac),
      average_0th_deriv_of_q_(average_0th_deriv_of_q) {
  if (avg_tscale_frac_ <= 0.0) {
    ERROR(
        "The specified avg_timescale_frac, the fraction of the damping "
        "timescale over which averaging is performed, ("
        << avg_tscale_frac_ << ") must be positive");
  }
}

template <size_t DerivOrder>
Averager<DerivOrder>::Averager(Averager&& rhs) noexcept
    : avg_tscale_frac_(std::move(rhs.avg_tscale_frac_)),
      average_0th_deriv_of_q_(std::move(rhs.average_0th_deriv_of_q_)),
      averaged_values_(std::move(rhs.averaged_values_)),
      nullopt_(std::move(rhs.nullopt_)),
      times_(std::move(rhs.times_)),
      raw_qs_(std::move(rhs.raw_qs_)),
      weight_k_(std::move(rhs.weight_k_)),
      tau_k_(std::move(rhs.tau_k_)) {}

template <size_t DerivOrder>
Averager<DerivOrder>& Averager<DerivOrder>::operator=(Averager&& rhs) noexcept {
  if (this != &rhs) {
    avg_tscale_frac_ = std::move(rhs.avg_tscale_frac_);
    average_0th_deriv_of_q_ = std::move(rhs.average_0th_deriv_of_q_);
    averaged_values_ = std::move(rhs.averaged_values_);
    nullopt_ = std::move(rhs.nullopt_);
    times_ = std::move(rhs.times_);
    raw_qs_ = std::move(rhs.raw_qs_);
    weight_k_ = std::move(rhs.weight_k_);
    tau_k_ = std::move(rhs.tau_k_);
  }
  return *this;
}

template <size_t DerivOrder>
const std::optional<std::array<DataVector, DerivOrder + 1>>&
Averager<DerivOrder>::operator()(const double time) const noexcept {
  std::ostringstream os;
  // os << "Inside averager operator.\n";
  // os << "  times_.size() = " << times_.size() << "\n"
  //   << "  times_[0] = " << times_[0] << "\n";
  // Parallel::printf(os.str());
  if (times_.size() > DerivOrder and time == times_[0]) {
    return averaged_values_;
  }
  return nullopt_;
}

template <size_t DerivOrder>
void Averager<DerivOrder>::clear() noexcept {
  averaged_values_ = std::nullopt;
  times_.clear();
  raw_qs_.clear();
  weight_k_ = 0.0;
  tau_k_ = 0.0;
}

template <size_t DerivOrder>
void Averager<DerivOrder>::update(const double time, const DataVector& raw_q,
                                  const DataVector& timescales) noexcept {
  if (not raw_qs_.empty()) {
    if (UNLIKELY(raw_q.size() != raw_qs_[0].size())) {
      ERROR("The number of components in the raw_q provided ("
            << raw_q.size()
            << ") does not match the size of previously supplied raw_q ("
            << raw_qs_[0].size() << ").");
    }
  } else {
    // This is the first call to update: initialize averaged values, weights and
    // effective time (with proper number of components)
    averaged_values_ =
        make_array<DerivOrder + 1>(DataVector(raw_q.size(), 0.0));
    weight_k_ = 0.0;
    tau_k_ = 0.0;
  }

  // Check if we need to update
  if (time == times_[0]) {
    return;
  }

  // Ensure that the number of timescales matches the number of components
  if (UNLIKELY(timescales.size() != raw_q.size())) {
    ERROR("The number of supplied timescales ("
          << timescales.size() << ") does not match the number of components ("
          << raw_q.size() << ").");
  }
  // Get the minimum damping time from all component timescales. This will be
  // used to determine the averaging timescale for ALL components.
  const double min_timescale = min(timescales);

  // Do not allow updates at or before last update time
  if (UNLIKELY(not times_.empty() and time <= last_time_updated())) {
    ERROR("The specified time t=" << time
                                  << " is at or before the last time "
                                     "updated, t_update="
                                  << last_time_updated() << ".");
  }

  // update deques
  times_.emplace_front(time);
  raw_qs_.emplace_front(raw_q);
  if (times_.size() > DerivOrder + 1) {
    // get rid of old data once we have a sufficient number of points
    times_.pop_back();
    raw_qs_.pop_back();
  }

  if (times_.size() > DerivOrder) {
    // we now have enough data points to begin averaging
    const double tau_avg = min_timescale * avg_tscale_frac_;
    const double tau_m = times_[0] - times_[1];

    // update the weights and effective time
    const double old_weight = weight_k_;
    weight_k_ = (tau_m + old_weight) * tau_avg / (tau_m + tau_avg);
    tau_k_ = (time * tau_m + tau_k_ * old_weight) * tau_avg /
             (weight_k_ * (tau_m + tau_avg));

    std::array<DataVector, DerivOrder + 1> raw_derivs = get_derivs();
    // use raw value if not using `average_0th_deriv_of_q`
    size_t start_ind = 0;
    if (not average_0th_deriv_of_q_) {
      (*averaged_values_)[0] = raw_qs_[0];
      start_ind = 1;
    }

    for (size_t i = start_ind; i <= DerivOrder; i++) {
      gsl::at(*averaged_values_, i) =
          (gsl::at(raw_derivs, i) * tau_m +
           gsl::at(*averaged_values_, i) * old_weight) *
          tau_avg / (weight_k_ * (tau_m + tau_avg));
    }
  }
}

template <size_t DerivOrder>
double Averager<DerivOrder>::last_time_updated() const noexcept {
  if (UNLIKELY(times_.empty())) {
    ERROR("The time history has not been updated yet.");
  }
  return times_[0];
}

template <size_t DerivOrder>
double Averager<DerivOrder>::average_time(const double time) const noexcept {
  if (LIKELY(times_.size() > DerivOrder and time == times_[0])) {
    return tau_k_;
  }
  ERROR(
      "Cannot return averaged values because the averager does not have "
      "up-to-date information.");
}

template <size_t DerivOrder>
std::array<DataVector, DerivOrder + 1> Averager<DerivOrder>::get_derivs()
    const noexcept {
  // initialize the finite difference coefs
  std::array<std::array<double, DerivOrder + 1>, DerivOrder + 1> coefs{};

  if constexpr (DerivOrder == 0) {
    // set coefs for function value
    coefs[0] = {{1.0}};
  } else if constexpr (DerivOrder == 1) {
    const double delta_t = times_[0] - times_[1];

    // set coefs for function value
    coefs[0] = {{1.0, 0.0}};
    // first deriv coefs
    coefs[1] = {{-1.0 / delta_t, -1.0 / delta_t}};
  } else if constexpr (DerivOrder == 2) {
    const double delta_t1 = times_[0] - times_[1];
    const double delta_t2 = times_[1] - times_[2];
    const double delta_t = delta_t1 + delta_t2;
    const double denom = 1.0 / delta_t;
    const double dts = delta_t1 * delta_t2;

    // set coefs for function value
    coefs[0] = {{1.0, 0.0, 0.0}};
    // first deriv coefs
    coefs[1] = {{(2.0 + delta_t2 / delta_t1) * denom, -delta_t / dts,
                 delta_t1 / delta_t2 * denom}};
    // second deriv coefs
    coefs[2] = {{2.0 / delta_t1 * denom, -2.0 / dts, 2.0 / delta_t2 * denom}};
  } else if constexpr (DerivOrder == 3) {
    const double delta_t1 = times_[1] - times_[0];
    const double delta_t2 = times_[2] - times_[0];
    const double delta_t3 = times_[3] - times_[0];
    const double dt1_minus_dt2 = times_[1] - times_[2];
    const double dt1_minus_dt3 = times_[1] - times_[3];
    const double dt2_minus_dt3 = times_[2] - times_[3];

    // set coefs for function value
    coefs[0] = {{1.0, 0.0, 0.0, 0.0}};
    // first deriv coefs
    coefs[1] = {
        {delta_t2 * delta_t3 / (delta_t1 * dt1_minus_dt2 * dt1_minus_dt3),
         -delta_t1 * delta_t3 / (delta_t2 * dt1_minus_dt2 * dt2_minus_dt3),
         delta_t1 * delta_t2 / (delta_t3 * dt1_minus_dt3 * dt2_minus_dt3),
         -1.0 / delta_t1 - 1.0 / delta_t2 - 1.0 / delta_t3}};
    // second deriv coefs
    coefs[2] = {{-2.0 * (delta_t2 + delta_t3) /
                     (delta_t1 * dt1_minus_dt2 * dt1_minus_dt3),
                 2.0 * (delta_t1 + delta_t3) /
                     (delta_t2 * dt1_minus_dt2 * dt2_minus_dt3),
                 -2.0 * (delta_t1 + delta_t2) /
                     (delta_t3 * dt1_minus_dt3 * dt2_minus_dt3),
                 2.0 * (delta_t1 + delta_t2 + delta_t3) /
                     (delta_t1 * delta_t2 * delta_t3)}};
    // third deriv coefs
    coefs[3] = {{6.0 / (delta_t1 * dt1_minus_dt2 * dt1_minus_dt3),
                 -6.0 / (delta_t2 * dt1_minus_dt2 * dt2_minus_dt3),
                 6.0 / (delta_t3 * dt1_minus_dt3 * dt2_minus_dt3),
                 -6.0 / (delta_t1 * delta_t2 * delta_t3)}};
  } else {
    ERROR("Cannot compute derivatives for a deriv order greater than 3.");
  }

  // compute derivatives
  std::array<DataVector, DerivOrder + 1> result =
      make_array<DerivOrder + 1>(DataVector(raw_qs_[0].size(), 0.0));
  for (size_t i = 0; i < DerivOrder + 1; i++) {
    for (size_t j = 0; j < DerivOrder + 1; j++) {
      gsl::at(result, i) += raw_qs_[j] * gsl::at(gsl::at(coefs, i), j);
    }
  }
  return result;
}

template <size_t DerivOrder>
void Averager<DerivOrder>::pup(PUP::er& p) noexcept {
  p | avg_tscale_frac_;
  p | average_0th_deriv_of_q_;
  p | averaged_values_;
  p | times_;
  p | raw_qs_;
  p | weight_k_;
  p | tau_k_;
}

template <size_t DerivOrder>
std::string Averager<DerivOrder>::get_output() const noexcept {
  std::ostringstream os;
  os << "avg_tscale_frac_ = " << avg_tscale_frac_ << "\n"
     << "average_0th_deriv_of_q_ = " << average_0th_deriv_of_q_ << "\n";
  if (averaged_values_.has_value()) {
    os << "Size of averaged_values_ = " << averaged_values_.value().size()
       << "\n";
  }
  os << "Size of times_ = " << times_.size() << "\n"
     << "times_ = (";
  for (size_t i = 0; i < times_.size() - 1; i++) {
    os << times_[i] << ", ";
  }
  os << times_[times_.size() - 1] << ")\n";
  os << "Size of raw_qs_ = " << raw_qs_.size() << "\n"
     << "weight_k_ = " << weight_k_ << "\n"
     << "tau_k_ = " << tau_k_;
  return os.str();
}

template <size_t DerivOrder>
bool operator==(const Averager<DerivOrder>& avg1,
                const Averager<DerivOrder>& avg2) {
  return (avg1.avg_tscale_frac_ == avg2.avg_tscale_frac_) and
         (avg1.average_0th_deriv_of_q_ == avg2.average_0th_deriv_of_q_) and
         (avg1.averaged_values_ == avg2.averaged_values_) and
         (avg1.times_ == avg2.times_) and (avg1.raw_qs_ == avg2.raw_qs_) and
         (avg1.weight_k_ == avg2.weight_k_) and (avg1.tau_k_ == avg2.tau_k_);
}

template <size_t DerivOrder>
bool operator!=(const Averager<DerivOrder>& avg1,
                const Averager<DerivOrder>& avg2) {
  return not(avg1 == avg2);
}

// explicit instantiations for deriv order (0,1,2,3)
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                            \
  template class Averager<DIM(data)>;                   \
  template bool operator==(const Averager<DIM(data)>&,  \
                           const Averager<DIM(data)>&); \
  template bool operator!=(const Averager<DIM(data)>&,  \
                           const Averager<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (0, 1, 2, 3))

#undef INSTANIATE
#undef DIM
