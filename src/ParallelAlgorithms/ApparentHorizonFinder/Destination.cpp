// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/ApparentHorizonFinder/Destination.hpp"

#include <string>

#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace ah {
std::ostream& operator<<(std::ostream& os, const Destination ordering) {
  switch (ordering) {
    case Destination::Observation:
      return os << "Observation";
    case Destination::ControlSystem:
      return os << "ControlSystem";
    default:
      ERROR("Unknown Destination type");
  }
}
}  // namespace ah

template <>
ah::Destination Options::create_from_yaml<ah::Destination>::create<void>(
    const Options::Option& options) {
  const auto ordering = options.parse_as<std::string>();
  if (ordering == "Observation") {
    return ah::Destination::Observation;
  } else if (ordering == "ControlSystem") {
    return ah::Destination::ControlSystem;
  }
  PARSE_ERROR(options.context(),
              "Destination must be 'Observation' or 'ControlSystem'");
}
