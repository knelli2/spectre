// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Parallel/Protocols/ArrayElementsAllocator.hpp"
#include "Utilities/ProtocolHelpers.hpp"

namespace intrp2::protocols {
/*!
 * \brief A protocol for compile time interpolation options in the metavariables
 * of an executable.
 *
 * A class conforming to this protocol is placed in the metavariables to choose
 * compile-time options for interpolation. The conforming class must:
 * - provide a type alias `elements_allocator` which conforms to the
 *   `Parallel::protocols::ArrayElementsAllocator` that will initialize the
 *   elements of the interpolation array component.
 * - provide a type alias `element_initializer` which is a simple action that
 *   takes no arguments which will initialize the `intrp2::Tags::DbAccess` tag
 *   of the DataBox. Note that the array index of this simple action must be a
 *   `std::string`.
 * - be named `intrp`.
 *
 * TODO: Add example
 */
struct Metavariables {
  template <typename ConformingType>
  struct test {
    using elements_allocator = typename ConformingType::elements_allocator;
    static_assert(
        tt::assert_conforms_to_v<elements_allocator,
                                 Parallel::protocols::ArrayElementsAllocator>);
    using element_initializer = typename ConformingType::element_initializer;
  };
};
}  // namespace intrp2::protocols
