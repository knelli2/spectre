// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system {
template <typename ControlSystem, typename ReturnTags, typename ArgumentTags>
struct ComputeControlError;

/*!
 * \brief DataBox mutator for calculating the control error of the given
 * `ControlSystem`
 *
 * \details Intended to be used in a `db::mutate_apply` call. Control errors may
 * have to edit other tags in the DataBox as they compute their control error
 * (like for size control). This acts as a safe way for the control error to
 * edit whatever tags it wants in the DataBox, while also being able to edit
 * itself. This is used as a safe alternative to something unsafe like
 *
 * \code {.cpp}
 * auto& control_error =
 *   db::get_mutable_reference<ControlError>(make_not_null(&box));
 * auto& tag1 = db::get_mutable_reference<Tag1>(make_not_null(&box));
 * auto& tag2 = db::get_mutable_reference<Tag2>(make_not_null(&box));
 *
 * control_error(make_not_null(&tag1), make_not_null(&tag1), args...);
 * \endcode
 *
 * \note The `ReturnTags` and `ArgumentTags` template parameters must be the
 * same as the `return_tags` and `argument_tags` type aliases of the control
 * error in the `ControlSystem`. They are passed in to the struct because they
 * are unable to be deduced in the `apply()` function itself.
 */
template <typename ControlSystem, typename... ReturnTags,
          typename... ArgumentTags>
struct ComputeControlError<ControlSystem, tmpl::list<ReturnTags...>,
                           tmpl::list<ArgumentTags...>> {
  using return_tags =
      tmpl::list<control_system::Tags::ControlError<ControlSystem>,
                 ReturnTags...>;
  using argument_tags = tmpl::list<ArgumentTags...>;

  static_assert(
      std::is_same_v<tmpl::list<ReturnTags...>,
                     typename ControlSystem::control_error::return_tags>,
      "The ReturnTags passed to ComputeControlError must be the same as the "
      "`return_tags` type alias of the ControlSystem::control_error struct.");
  static_assert(
      std::is_same_v<tmpl::list<ArgumentTags...>,
                     typename ControlSystem::control_error::argument_tags>,
      "The ArgumentTags passed to ComputeControlError must be the same as the "
      "`argument_tags` type alias of the ControlSystem::control_error struct.");

  template <typename Metavariables, typename... TupleTags>
  static DataVector apply(
      const gsl::not_null<typename ControlSystem::control_error*> control_error,
      const gsl::not_null<typename ReturnTags::type*>... return_tags,
      const typename ArgumentTags::type&... argument_tags,
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& function_of_time_name,
      const tuples::TaggedTuple<TupleTags...>& measurements) {
    return control_error->operator()(return_tags..., argument_tags..., cache,
                                     time, function_of_time_name, measurements);
  }
};
}  // namespace control_system
