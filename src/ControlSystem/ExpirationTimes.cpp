// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ExpirationTimes.hpp"

#include "DataStructures/DataVector.hpp"

namespace control_system {
double function_of_time_expiration_time(
    const double time, const DataVector& measurement_timescales,
    const int measurements_per_update) {
  return time + measurements_per_update * min(measurement_timescales);
}

double measurement_expiration_time(const double time,
                                   const DataVector& measurement_timescales,
                                   const int measurements_per_update) {
  return function_of_time_expiration_time(time, measurement_timescales,
                                          measurements_per_update) -
         0.5 * min(measurement_timescales);
}
}  // namespace control_system
