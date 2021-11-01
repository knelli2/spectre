// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "ControlSystem/Component.hpp"
#include "ControlSystem/Tags.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/OptionTags.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Time/Tags.hpp"
#include "Utilities/TypeTraits.hpp"

namespace {

template <bool Override>
struct MetavariablesWithoutControlSystem {
  static constexpr size_t volume_dim = 3;
  static constexpr bool override_functions_of_time = Override;
  using component_list = tmpl::list<>;
};

using one_system = control_system::TestHelpers::System<
    2, control_system::TestHelpers::TestStructs_detail::LabelA,
    control_system::TestHelpers::Measurement<
        control_system::TestHelpers::TestStructs_detail::LabelA>>;

template <bool Override>
struct MetavariablesWithControlSystem {
  static constexpr size_t volume_dim = 3;
  static constexpr bool override_functions_of_time = Override;
  using control_systems = tmpl::list<one_system>;
  using control_components =
      control_system::control_components<MetavariablesWithControlSystem,
                                         control_systems>;

  using component_list = tmpl::flatten<tmpl::list<control_components>>;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.FunctionsOfTime.Tags", "[Domain][Unit]") {
  TestHelpers::db::test_simple_tag<domain::Tags::FunctionsOfTime>(
      "FunctionsOfTime");

  CHECK(std::is_same_v<
        domain::Tags::FunctionsOfTime::option_tags<
            MetavariablesWithoutControlSystem<true>>,
        tmpl::list<
            domain::OptionTags::DomainCreator<
                MetavariablesWithoutControlSystem<true>::volume_dim>,
            domain::FunctionsOfTime::OptionTags::FunctionOfTimeFile,
            domain::FunctionsOfTime::OptionTags::FunctionOfTimeNameMap>>);

  CHECK(std::is_same_v<
        domain::Tags::FunctionsOfTime::option_tags<
            MetavariablesWithoutControlSystem<false>>,
        tmpl::list<domain::OptionTags::DomainCreator<
                       MetavariablesWithoutControlSystem<false>::volume_dim>,
                   OptionTags::InitialTimeStep>>);

  CHECK(std::is_same_v<
        domain::Tags::FunctionsOfTime::option_tags<
            MetavariablesWithControlSystem<false>>,
        tmpl::list<
            domain::OptionTags::DomainCreator<
                MetavariablesWithControlSystem<false>::volume_dim>,
            OptionTags::InitialTimeStep,
            control_system::OptionTags::ControlSystemInputs<one_system>>>);
}
