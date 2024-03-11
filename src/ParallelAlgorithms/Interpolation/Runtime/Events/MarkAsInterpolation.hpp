// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>

#include "DataStructures/DataBox/Access.hpp"

namespace intrp2::Events {
/// To designate an Event as an interpolation event. Necessary to be able to
/// initialize the parallel component properly
struct MarkAsInterpolation {
  // This will be used to initialize the Accesses that correspond to the
  // individual targets
  virtual std::unique_ptr<db::Access> initialize_target_element_box() const = 0;
  virtual ~MarkAsInterpolation() {}
};
}  // namespace intrp2::Events
