// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "ControlSystem/Protocols/Measurement.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace control_system {
struct TrackTranslation : tt::ConformsTo<protocols::Measurement> {
  using submeasurements = tmpl::list<SubTrackTranslation>;
};
}  // namespace control_system
