// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Targets/AngularOrdering.hpp"

#include <string>

#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace intrp {
std::ostream& operator<<(std::ostream& os, const AngularOrdering ordering) {
  switch (ordering) {
    case AngularOrdering::Strahlkorper:
      return os << "Strahlkorper";
    case AngularOrdering::CCE:
      return os << "CCE";
    default:
      ERROR("Unknown AngularOrdering type");
  }
}
}  // namespace intrp

template <>
intrp::AngularOrdering
Options::create_from_yaml<intrp::AngularOrdering>::create<void>(
    const Options::Option& options) {
  const auto ordering = options.parse_as<std::string>();
  if (ordering == "Strahlkorper") {
    return intrp::AngularOrdering::Strahlkorper;
  } else if (ordering == "CCE") {
    return intrp::AngularOrdering::CCE;
  }
  PARSE_ERROR(options.context(),
              "AngularOrdering must be 'Strahlkorper' or 'CCE'");
}
