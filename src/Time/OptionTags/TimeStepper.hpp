// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <string>

#include "Evolution/Tags.hpp"
#include "Options/String.hpp"

namespace OptionTags {
/// \ingroup OptionTagsGroup
/// \ingroup TimeGroup
template <typename StepperType>
struct TimeStepper {
  static std::string name() { return "TimeStepper"; }
  inline const static std::string help{"The time stepper"};
  using type = std::unique_ptr<StepperType>;
  using group = evolution::OptionTags::Group;
};
}  // namespace OptionTags
