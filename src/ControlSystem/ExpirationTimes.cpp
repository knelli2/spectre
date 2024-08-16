// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ExpirationTimes.hpp"

#include "DataStructures/DataVector.hpp"

namespace control_system {
double function_of_time_expiration_time(
    const double time, const double fraction,
    const DataVector& old_measurement_timescales,
    const DataVector& new_measurement_timescales,
    const int measurements_per_update) {
  return time +
         fraction * (min(old_measurement_timescales) +
                     min(new_measurement_timescales)) +
         (measurements_per_update - 1) * min(new_measurement_timescales);
}

double measurement_expiration_time(const double time, const double fraction,
                                   const DataVector& old_measurement_timescales,
                                   const DataVector& new_measurement_timescales,
                                   const int measurements_per_update) {
  return function_of_time_expiration_time(
             time, fraction, old_measurement_timescales,
             new_measurement_timescales, measurements_per_update) -
         fraction * 0.5 * min(new_measurement_timescales);
}
}  // namespace control_system
