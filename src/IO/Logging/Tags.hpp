// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"
#include "Utilities/PrettyType.hpp"

/// \cond
enum class Verbosity;
/// \endcond

/// \ingroup LoggingGroup
/// Items related to logging
namespace logging {
namespace OptionTags {
/// \ingroup OptionTagsGroup
/// \ingroup LoggingGroup
template <typename OptionsGroup>
struct Verbosity {
  using type = ::Verbosity;
  inline const static std::string help{"Verbosity"};
  using group = OptionsGroup;
};
}  // namespace OptionTags

namespace Tags {
/// \ingroup LoggingGroup
/// \brief Tag for putting `::Verbosity` in a DataBox.
template <typename OptionsGroup>
struct Verbosity : db::SimpleTag {
  using type = ::Verbosity;
  static std::string name() {
    return "Verbosity(" + pretty_type::name<OptionsGroup>() + ")";
  }

  using option_tags = tmpl::list<OptionTags::Verbosity<OptionsGroup>>;
  static constexpr bool pass_metavariables = false;
  static ::Verbosity create_from_options(const ::Verbosity& verbosity) {
    return verbosity;
  }
};
}  // namespace Tags
}  // namespace logging
