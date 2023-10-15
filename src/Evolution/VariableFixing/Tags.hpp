// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"
#include "Utilities/PrettyType.hpp"

namespace OptionTags {
/*!
 * \ingroup OptionGroupsGroup
 * \brief Groups the variable fixer configurations in the input file.
 */
struct VariableFixingGroup {
  static std::string name() { return "VariableFixing"; }
  inline const static std::string help {"Options for variable fixing"};
};

/*!
 * \ingroup OptionTagsGroup
 * \brief The global cache tag that retrieves the parameters for the variable
 * fixer from the input file.
 */
template <typename VariableFixerType>
struct VariableFixer {
  inline const static std::string help {"Options for the variable fixer"};
  using type = VariableFixerType;
  static std::string name() { return pretty_type::name<VariableFixerType>(); }
  using group = VariableFixingGroup;
};
}  // namespace OptionTags

namespace Tags {
/*!
 * \brief The global cache tag for the variable fixer
 */
template <typename VariableFixerType>
struct VariableFixer : db::SimpleTag {
  using type = VariableFixerType;
  using option_tags =
      tmpl::list<::OptionTags::VariableFixer<VariableFixerType>>;

  static constexpr bool pass_metavariables = false;
  static VariableFixerType create_from_options(
      const VariableFixerType& variable_fixer) {
    return variable_fixer;
  }
};
}  // namespace Tags
