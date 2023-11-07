// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Falloff.hpp"

#include <ostream>
#include <string>

#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {

std::ostream& operator<<(std::ostream& os, const Falloff falloff) {
  switch (falloff) {
    case Falloff::Linear:
      return os << "Linear";
    case Falloff::Inverse:
      return os << "Inverse";
    default:
      ERROR(
          "Unknown "
          "domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff type");
  }
}
}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions

template <>
domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff
Options::create_from_yaml<
    domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff>::
    create<void>(const Options::Option& options) {
  const auto falloff = options.parse_as<std::string>();
  if (falloff == "Linear") {
    return domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff::Linear;
  } else if (falloff == "Inverse") {
    return domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff::
        Inverse;
  }
  PARSE_ERROR(options.context(), "Falloff must be 'Linear' or 'Inverse'");
}
