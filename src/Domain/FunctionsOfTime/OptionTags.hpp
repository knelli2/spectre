// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>
#include <unordered_map>

#include "Options/Auto.hpp"
#include "Options/String.hpp"

namespace domain::FunctionsOfTime::OptionTags {
/*!
 * \ingroup OptionGroupsGroup
 * \brief Groups options for reading in FunctionOfTime data from SpEC
 */
struct CubicFunctionOfTimeOverride {
  inline const static std::string help{
      "Options for importing FunctionOfTimes from SpEC"};
};

/*!
 * \brief Path to an H5 file containing SpEC FunctionOfTime data
 */
struct FunctionOfTimeFile {
  using type = Options::Auto<std::string, Options::AutoLabel::None>;
  inline const static std::string help{
      "Path to an H5 file containing SpEC FunctionOfTime data"};
  using group = CubicFunctionOfTimeOverride;
};

/*!
 * \brief Pairs of strings mapping SpEC FunctionOfTime names to SpECTRE names
 */
struct FunctionOfTimeNameMap {
  using type = std::map<std::string, std::string>;
  inline const static std::string help{
      "String pairs mapping spec names to spectre names"};
  using group = CubicFunctionOfTimeOverride;
};
}  // namespace domain::FunctionsOfTime::OptionTags
