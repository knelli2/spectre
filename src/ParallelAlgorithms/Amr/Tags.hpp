// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "Options/String.hpp"

/// Options for AMR
namespace amr::OptionTags {

struct AmrGroup {
  static std::string name() { return "Amr"; }
  inline const static std::string help
      {"Options for adaptive mesh refinement (AMR)"};
};

}  // namespace amr::OptionTags
