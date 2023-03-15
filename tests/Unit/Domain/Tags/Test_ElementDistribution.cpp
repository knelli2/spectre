// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "Framework/TestCreation.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Domain/ElementDistribution.hpp"
#include "Domain/Tags/ElementDistribution.hpp"

SPECTRE_TEST_CASE("Unit.Domain.Tags.ElementDistribution",
                  "[Unit][Domain]") {
  TestHelpers::db::test_simple_tag<domain::Tags::ElementWeight>("ElementWeight");
  CHECK(TestHelpers::test_option_tag<domain::OptionTags::ElementWeight>(
            "Uniform") == domain::ElementWeight::Uniform);
}
