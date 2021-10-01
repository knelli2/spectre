// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "DataStructures/DataVector.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system {
namespace ControlErrors {
struct Translation {
  template <size_t DerivOrder, typename... TupleTags>
  static DataVector apply(
      const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& integral_vector =
        get<control_system::QueueTags::MeasureTranslation>(measurements);
    DataVector target_value{DerivOrder - 1, 0.0};
    // integral(psi^2 * x_i) / integral(psi^2)
    // basically a CoM
    for (size_t i = 1; i < DerivOrder; i++) {
      target_value[i - 1] = integral_vector[i] / integral_vector[0];
    }

    const double& center_of_interval = integral_vector[DerivOrder];
    return target_value - center_of_interval;
  }
};
}  // namespace ControlErrors
}  // namespace control_system
