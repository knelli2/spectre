// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Falloff.hpp"
#include "Framework/TestCreation.hpp"
#include "Utilities/GetOutput.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {

SPECTRE_TEST_CASE("Unit.Domain.CoordinateMaps.Shape.Falloff",
                  "[Domain][Unit]") {
  CHECK(get_output(Falloff::Linear) == "Linear");
  CHECK(get_output(Falloff::Inverse) == "Inverse");
  CHECK(TestHelpers::test_creation<Falloff>("Linear") == Falloff::Linear);
  CHECK(TestHelpers::test_creation<Falloff>("Inverse") == Falloff::Inverse);
}

}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
