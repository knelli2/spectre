// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ExpirationTimes.hpp"

#include "DataStructures/DataVector.hpp"

namespace control_system {
double function_of_time_expiration_time(
    const double time, const double fraction,
    const DataVector& old_measurement_timescales,
    const DataVector& new_measurement_timescales,
    const int measurements_per_update,
    const ExpirationMethods expiration_method, const bool initial) {
  if (expiration_method == ExpirationMethods::spec) {
    // If we are trying to calculate the initial expiration time, we set the
    // expiration time earlier because the first measurement counts towards the
    // initial set of measurements whereas in subsequent updates, this "first"
    // measurement when we compute the update and expiration times doesn't count
    // towards the current set of measurements
    return time +
           static_cast<double>(measurements_per_update - (initial ? 1 : 0)) *
               min(new_measurement_timescales);
  } else {
    return time + min(old_measurement_timescales) +
           (static_cast<double>(measurements_per_update - 1 -
                                (initial ? 1 : 0)) +
            fraction) *
               min(new_measurement_timescales);
  }
}

double measurement_expiration_time(const double time, const double fraction,
                                   const DataVector& old_measurement_timescales,
                                   const DataVector& new_measurement_timescales,
                                   const int measurements_per_update,
                                   const ExpirationMethods expiration_method,
                                   const bool initial) {
  if (expiration_method == ExpirationMethods::spec) {
    return function_of_time_expiration_time(
               time, fraction, old_measurement_timescales,
               new_measurement_timescales, measurements_per_update,
               expiration_method, initial) -
           0.5 * min(new_measurement_timescales);
  } else {
    return function_of_time_expiration_time(
               time, fraction, old_measurement_timescales,
               new_measurement_timescales, measurements_per_update,
               expiration_method, initial) -
           (fraction - 0.5) * min(new_measurement_timescales);
  }
}
}  // namespace control_system
