// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system::protocols {
/*!
 * \brief Definition of a control error
 *
 * A control error is used within a control system to compute how far off the
 * the value you are controlling is from its expected value.
 *
 * A conforming type must specify:
 *
 * - a `static constexpr size_t expected_number_of_excisions` which specifies
 *   the number of excisions necessary in order to compute the control error.
 * - a type alias `object_centers` to a `domain::object_list` of
 *   `domain::ObjectLabel`s. These are the objects that will require the
 *   `domain::Tags::ObjectCenter`s tags to be in the GlobalCache for this
 *   control system to work.
 * - a `return_tags` type alias to a `tmpl::list` of DataBox tags to be passed
 *   to the `operator()` so the ControlError can be used in a
 *   `db::mutate_apply`
 * - an `argument_tags` type alias to a `tmpl::list` of DataBox tags to be
 *    passed to the `operator()` so the ControlError can be used in a
 *   `db::mutate_apply`
 * - a call operator that returns a DataVector. This call operator show have
 *   the following signature
 *
 * \code {.cpp}
 * DataVector operator()(
 *     const gsl::not_null<ReturnTags::type...*> return_tags,
 *     const ArgumentTags::type&... argument_Tags,
 *     const Parallel::GlobalCache<Metavariables>& cache,
 *     const double time,
 *     const std::string& function_of_time_name,
 *     const tuples::TaggedTuple<TupleArgs...>& measurements)
 * \endcode
 *
 * where the `ReturnTags` and `ArgumentTags` correspond to the type aliases. If
 * there are no return tags, then they can be left out of the function
 * signature. Same for the argument tags. However, the `cache`, `time`,
 * `function_of_time_name`, and `measurements` must always be present. Here is
 * an example:
 *
 *   \snippet Helpers/ControlSystem/Examples.hpp ControlError
 */
struct ControlError {
  template <typename ConformingType>
  struct test {
    struct DummyMetavariables;
    struct DummyTupleTags;

    static constexpr size_t expected_number_of_excisions =
        ConformingType::expected_number_of_excisions;

    using object_centers = typename ConformingType::object_centers;

    using return_tags = typename ConformingType::return_tags;
    using argument_tags = typename ConformingType::argument_tags;
  };
};
}  // namespace control_system::protocols
