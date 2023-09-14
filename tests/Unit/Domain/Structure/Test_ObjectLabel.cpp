// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Domain/Structure/ObjectLabel.hpp"
#include "Utilities/GetOutput.hpp"

SPECTRE_TEST_CASE("Unit.Domain.ObjectLabel", "[Domain][Unit]") {
  CHECK(name(domain::ObjectLabel::A) == "A");
  CHECK(get_output(domain::ObjectLabel::A) == "A");
  CHECK(name(domain::ObjectLabel::B) == "B");
  CHECK(get_output(domain::ObjectLabel::B) == "B");
  CHECK(name(domain::ObjectLabel::C) == "C");
  CHECK(get_output(domain::ObjectLabel::C) == "C");
  CHECK(name(domain::ObjectLabel::None) == "");
  CHECK(get_output(domain::ObjectLabel::None) == "");

  // This is a regression test to ensure nobody accidentally changes the values
  // of the enums
  CHECK(static_cast<int>(domain::ObjectLabel::A) == 0);
  CHECK(static_cast<int>(domain::ObjectLabel::B) == 1);
  CHECK(static_cast<int>(domain::ObjectLabel::C) == 2);
  CHECK(static_cast<int>(domain::ObjectLabel::None) == -1);
}
