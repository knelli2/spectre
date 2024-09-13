// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "Framework/TestCreation.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "Utilities/GetOutput.hpp"

namespace intrp2 {
SPECTRE_TEST_CASE(
    "Unit.NumericalAlgorithms.InterpolationTarget.Runtime.AngularOrdering",
    "[Unit]") {
  CHECK(get_output(AngularOrdering::Cce) == "Cce");
  CHECK(get_output(AngularOrdering::Strahlkorper) == "Strahlkorper");
  CHECK(TestHelpers::test_creation<AngularOrdering>("Cce") ==
        AngularOrdering::Cce);
  CHECK(TestHelpers::test_creation<AngularOrdering>("Strahlkorper") ==
        AngularOrdering::Strahlkorper);
}
}  // namespace intrp2
