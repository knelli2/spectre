// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system::protocols {
/// \brief Definition of a control error
///
/// A control error is used within a control system to compute how far off the
/// the value you are controlling is away from its expected value.
///
/// A conforming type must specify:
///
/// - an operator that returns a DataVector with a signature the same as in the
///   example shown here:
///   \snippet Helpers/ControlSystem/Examples ControlError
struct ControlError {
  template <typename ConformingType>
  struct test {
    struct DummyMetavariables;
    struct DummyTupleTags;

    static_assert(
        std::is_same_v<
            DataVector,
            decltype(ConformingType{}(
                std::declval<
                    const Parallel::GlobalCache<DummyMetavariables>&>(),
                std::declval<const double>(),
                std::declval<const std::string&>(),
                std::declval<const tuples::TaggedTuple<DummyTupleTags>&>()))>);
  };
};
}  // namespace control_system::protocols
