// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

namespace intrp2::protocols {
/*!
 * \brief A protocol for compile time interpolation options in the metavariables
 * of an executable.
 *
 * A class conforming to this protocol is placed in the metavariables to choose
 * compile-time options for interpolation. The conforming class must:
 * - provide a static member variable `bool include_dense_triggers`: Whether or
 *   not the executable uses dense triggers or not. This is necessary so the
 *   interpolation parallel component can be constructed from dense triggers as
 *   well as standard triggers.
 * - be names `intrp`.
 *
 * TODO: Add example
 */
struct Metavariables {
  template <typename ConformingType>
  struct test {
    using include_dense_triggers_type = const bool;
    using include_dense_triggers_return_type =
        decltype(ConformingType::include_dense_triggers);
    static_assert(std::is_same_v<include_dense_triggers_type,
                                 include_dense_triggers_return_type>,
                  "The metavariable 'enable_time_dependent_maps' should be a "
                  "static constexpr bool.");
  };
};
}  // namespace intrp2::protocols
