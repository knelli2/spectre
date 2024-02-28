// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"

#include <string>

#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace intrp2 {
std::ostream& operator<<(std::ostream& os, const AngularOrdering ordering) {
  switch (ordering) {
    case AngularOrdering::Strahlkorper:
      return os << "Strahlkorper";
    case AngularOrdering::Cce:
      return os << "Cce";
    default:
      ERROR("Unknown AngularOrdering type");
  }
}
}  // namespace intrp2

template <>
intrp2::AngularOrdering
Options::create_from_yaml<intrp2::AngularOrdering>::create<void>(
    const Options::Option& options) {
  const auto ordering = options.parse_as<std::string>();
  if (ordering == "Strahlkorper") {
    return intrp2::AngularOrdering::Strahlkorper;
  } else if (ordering == "Cce") {
    return intrp2::AngularOrdering::Cce;
  }
  PARSE_ERROR(options.context(),
              "AngularOrdering must be 'Strahlkorper' or 'Cce'");
}
