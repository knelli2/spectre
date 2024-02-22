// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ControlErrors/Translation.hpp"

#include <cstddef>

namespace control_system::ControlErrors {
template struct Translation<1_st>;
template struct Translation<2_st>;
}  // namespace control_system::ControlErrors
