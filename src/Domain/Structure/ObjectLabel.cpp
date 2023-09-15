// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Structure/ObjectLabel.hpp"

#include <ostream>

#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Literals.hpp"

namespace domain {
std::string name(const ObjectLabel x) {
  if (x == ObjectLabel::A) {
    return "A"s;
  } else if (x == ObjectLabel::B) {
    return "B"s;
  } else if (x == ObjectLabel::C) {
    return "C"s;
  } else if (x == ObjectLabel::None) {
    return ""s;
  } else {
    ERROR("Unknown object label!");
  }
}

std::ostream& operator<<(std::ostream& s, const ObjectLabel x) {
  return s << name(x);
}
}  // namespace domain

template <>
domain::ObjectLabel
Options::create_from_yaml<domain::ObjectLabel>::create<void>(
    const Options::Option& options) {
  const auto ordering = options.parse_as<std::string>();
  if (ordering == "A") {
    return domain::ObjectLabel::A;
  } else if (ordering == "B") {
    return domain::ObjectLabel::B;
  } else if (ordering == "C") {
    return domain::ObjectLabel::C;
  } else if (ordering == "None") {
    return domain::ObjectLabel::None;
  }

  PARSE_ERROR(options.context(),
              "ObjectLabel must be 'A', 'B', 'C', or 'None'");
}
