// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp2::protocols {
/*!
 * \brief A protocol for the target of an interpolation.
 *
 * \details A struct conforming to the `Target` protocol must have
 *
 * - a type alias `temporal_id_tag` which is a tag in the element DataBox whose
 *   type will be used to distinguish interpolations at different times.
 *
 * - a type alias `frame` which is either the frame that the points are in, or
 *   `NoSuchType` to represent that the frame can be chosen at runtime.
 *
 * - a type alias `points` which is a struct that conforms to the
 *   `intrp2::protocols::Points` protocol.
 *
 * - a type alias `compile_time_callbacks` which is a list of structs that
 *   conform to the `intrp2::protocols::Callback` protocol. These callbacks will
 *   always be run after an interpolation and require no input from the user.
 *
 * - a type alias `possible_runtime_callbacks` which is a list of structs that
 *   conform to the `intrp2::protocols;:Callback` protocol. These callbacks will
 *   be constructed from the input file and the user can choose which to run at
 *   runtime.
 *
 * TODO: Add example
 */
struct Target {
  template <typename ConformingType>
  struct test {
    using temporal_id_tag = typename ConformingType::temporal_id_tag;
    using frame = typename ConformingType::frame;
    using points = typename ConformingType::points;
    static_assert(tt::assert_conforms_to_v<points, Points>);

    using possible_runtime_callbacks =
        typename ConformingType::possible_runtime_callbacks;
    static_assert(tmpl::all<possible_runtime_callbacks,
                            tt::assert_conforms_to<tmpl::_1, Callback>>::value);
    using compile_time_callbacks =
        typename ConformingType::compile_time_callbacks;
    static_assert(tmpl::all<compile_time_callbacks,
                            tt::assert_conforms_to<tmpl::_1, Callback>>::value);
  };
};
}  // namespace intrp2::protocols
