// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ExpirationTimes.hpp"

#include "DataStructures/DataVector.hpp"

namespace control_system {
double function_of_time_expiration_time(
    const double time, const DataVector& old_measurement_timescales,
    const DataVector& new_measurement_timescales,
    const int measurements_per_update) {
  double expiration_time = time + min(old_measurement_timescales);
  const double min_new_measurement_timescale = min(new_measurement_timescales);
  // We calculate the new expiration time in this way so it matches with how the
  // control system trigger calculates the next measurement time so doubles line
  // up and aren't different by roundoff.
  for (int i = 0; i < measurements_per_update; i++) {
    expiration_time += min_new_measurement_timescale;
  }
  return expiration_time;
}

double measurement_expiration_time(const double time,
                                   const DataVector& old_measurement_timescales,
                                   const DataVector& new_measurement_timescales,
                                   const int measurements_per_update) {
  return function_of_time_expiration_time(time, old_measurement_timescales,
                                          new_measurement_timescales,
                                          measurements_per_update) -
         0.5 * min(new_measurement_timescales);
}
}  // namespace control_system
