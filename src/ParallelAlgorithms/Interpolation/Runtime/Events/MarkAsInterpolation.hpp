// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>

#include "ParallelAlgorithms/Interpolation/Runtime/AccessWrapper.hpp"
#include "Utilities/Gsl.hpp"

namespace intrp2::Events {
/// To designate an Event as an interpolation event. Necessary to be able to
/// initialize the parallel component properly
struct MarkAsInterpolation {
  // This will be used to initialize the Accesses that correspond to the
  // individual targets
  // TODO: This needs to take an intrp2::AccessWrapper because it needs to be
  // created in the event because that knows the box type
  virtual void initialize_target_element_box(
      const gsl::not_null<std::unique_ptr<intrp2::AccessWrapper>*>
          access_wrapper) const = 0;
  virtual ~MarkAsInterpolation() {}
};
}  // namespace intrp2::Events
