// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ControlErrors/Translation.hpp"

#include <optional>
#include <pup.h>
#include <pup_stl.h>
#include <string>

#include "Options/Context.hpp"
#include "Options/ParseError.hpp"

namespace control_system::ControlErrors {
RadialTranslation::RadialTranslation(
    const std::array<double, 2>& radii,
    const std::optional<std::array<double, 2>>& safety_distances,
    const Options::Context& context)
    : inner_outer_radius_(radii),
      averaged_radius_(0.5 * (inner_outer_radius_[0] + inner_outer_radius_[1])),
      safety_distances_(safety_distances) {
  if (inner_outer_radius_[0] >= inner_outer_radius_[1]) {
    PARSE_ERROR(context, "Inner radius "
                             << inner_outer_radius_[0]
                             << " must be less than the outer radius "
                             << inner_outer_radius_[1]);
  }
}

void RadialTranslation::pup(PUP::er& p) {
  p | inner_outer_radius_;
  p | averaged_radius_;
  p | safety_distances_;
}
}  // namespace control_system::ControlErrors
