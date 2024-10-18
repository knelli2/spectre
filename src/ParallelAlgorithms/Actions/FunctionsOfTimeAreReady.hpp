// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/ArrayCollection/IsDgElementCollection.hpp"
#include "Parallel/ArrayCollection/PerformAlgorithmOnElement.hpp"
#include "Parallel/ArrayCollection/Tags/ElementLocations.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/Callback.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/Actions/GetItemFromDistributedObject.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/StdHelpers.hpp"

/// \cond
namespace Tags {
struct Time;
}  // namespace Tags
namespace domain::Tags {
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace tuples {
template <class... Tags>
class TaggedTuple;
}  // namespace tuples
namespace control_system::Tags {
struct Verbosity;
}  // namespace control_system::Tags
/// \endcond

namespace domain {
namespace detail {
template <typename CacheTag, typename Callback, typename Metavariables,
          typename ArrayIndex, typename Component, typename... Args>
bool functions_of_time_are_ready_impl(
    Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& array_index,
    const Component* /*meta*/, const double time,
    const std::optional<std::unordered_set<std::string>>& functions_to_check,
    Args&&... args) {
  if constexpr (Parallel::is_in_mutable_global_cache<Metavariables, CacheTag>) {
    const auto& proxy =
        ::Parallel::get_parallel_component<Component>(cache)[array_index];
    const Parallel::ArrayComponentId array_component_id =
        [&]() -> Parallel::ArrayComponentId {
      if constexpr (Parallel::is_dg_element_collection_v<Component>) {
        static_assert(
            std::is_same_v<std::tuple<ElementId<3>, double>,
                           std::tuple<std::decay_t<Args>...>>,
            "This currently assumes the only argument is the ElementId. If you "
            "need extra arguments, this can be generalized. We need the "
            "ElementId instead of the array_index (which is the nodegroup ID) "
            "since each ArrayComponentId is only allowed to register one "
            "callback, additional ones are ignored. This means of a "
            "DgElementCollection only 1 callback _per node_ would be "
            "registered, while we need 1 callback for each element. An "
            "alternative approach would be to have the callback be a broadcast "
            "to all elements instead of one specific one.");
        return Parallel::make_array_component_id<Component>(
            std::get<0>(std::tie(args...)));
      } else {
        return Parallel::make_array_component_id<Component>(array_index);
      }
    }();

    return Parallel::mutable_cache_item_is_ready<CacheTag>(
        cache, array_component_id,
        [&functions_to_check, &proxy, &time,
         &args...](const std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
                       functions_of_time) {
          using ::operator<<;
          ASSERT(
              alg::all_of(
                  functions_to_check.value_or(
                      std::unordered_set<std::string>{}),
                  [&functions_of_time](const std::string& function_to_check) {
                    return functions_of_time.count(function_to_check) == 1;
                  }),
              "Not all functions to check ("
                  << functions_to_check.value() << ") are in the global cache ("
                  << keys_of(functions_of_time) << ")");
          for (const auto& [name, f_of_t] : functions_of_time) {
            if (functions_to_check.has_value() and
                functions_to_check->count(name) == 0) {
              continue;
            }
            const double expiration_time = f_of_t->time_bounds()[1];
            if (time > expiration_time) {
              return std::unique_ptr<Parallel::Callback>(
                  new Callback(proxy, std::forward<Args>(args)...));
            }
          }
          return std::unique_ptr<Parallel::Callback>{};
        });
  } else {
    (void)cache;
    (void)array_index;
    (void)time;
    (void)functions_to_check;
    EXPAND_PACK_LEFT_TO_RIGHT((void)args);
    return true;
  }
}
}  // namespace detail

/// \ingroup ComputationalDomainGroup
/// Check that functions of time are up-to-date.
///
/// Check that functions of time in \p CacheTag with names in \p
/// functions_to_check are ready at time \p time.  If  \p functions_to_check is
/// a `std::nullopt`, checks all functions in \p CacheTag.  If any function is
/// not ready, schedules a `Parallel::PerformAlgorithmCallback` with the
/// GlobalCache..
template <typename CacheTag, typename Metavariables, typename ArrayIndex,
          typename Component>
bool functions_of_time_are_ready_algorithm_callback(
    Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& array_index,
    const Component* component_p, const double time,
    const std::optional<std::unordered_set<std::string>>& functions_to_check =
        std::nullopt) {
  using ProxyType =
      std::decay_t<decltype(::Parallel::get_parallel_component<Component>(
          cache)[array_index])>;
  return detail::functions_of_time_are_ready_impl<
      CacheTag, Parallel::PerformAlgorithmCallback<ProxyType>>(
      cache, array_index, component_p, time, functions_to_check);
}

/// \ingroup ComputationalDomainGroup
/// Check that functions of time are up-to-date.
///
/// Check that functions of time in \p CacheTag with names in \p
/// functions_to_check are ready at time \p time.  If  \p functions_to_check is
/// a `std::nullopt`, checks all functions in \p CacheTag.  If any function is
/// not ready, schedules a `Parallel::SimpleActionCallback` with the GlobalCache
/// which calls the simple action passed in as a template parameter. The `Args`
/// are forwareded to the callback.
template <typename CacheTag, typename SimpleAction, typename Metavariables,
          typename ArrayIndex, typename Component, typename... Args>
bool functions_of_time_are_ready_simple_action_callback(
    Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& array_index,
    const Component* component_p, const double time,
    const std::optional<std::unordered_set<std::string>>& functions_to_check,
    Args&&... args) {
  using ProxyType =
      std::decay_t<decltype(::Parallel::get_parallel_component<Component>(
          cache)[array_index])>;
  return detail::functions_of_time_are_ready_impl<
      CacheTag,
      Parallel::SimpleActionCallback<SimpleAction, ProxyType, Args...>>(
      cache, array_index, component_p, time, functions_to_check,
      std::forward<Args>(args)...);
}

/// \ingroup ComputationalDomainGroup
/// Check that functions of time are up-to-date.
///
/// Check that functions of time in \p CacheTag with names in \p
/// functions_to_check are ready at time \p time.  If  \p functions_to_check is
/// a `std::nullopt`, checks all functions in \p CacheTag.  If any function is
/// not ready, schedules a `Parallel::ThreadedActionCallback` with the
/// GlobalCache which calls the threaded action passed in as a template
/// parameter. The `Args` are forwareded to the callback.
template <typename CacheTag, typename ThreadedAction, typename Metavariables,
          typename ArrayIndex, typename Component, typename... Args>
bool functions_of_time_are_ready_threaded_action_callback(
    Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& array_index,
    const Component* component_p, const double time,
    const std::optional<std::unordered_set<std::string>>& functions_to_check,
    Args&&... args) {
  using ProxyType =
      std::decay_t<decltype(::Parallel::get_parallel_component<Component>(
          cache)[array_index])>;
  return detail::functions_of_time_are_ready_impl<
      CacheTag,
      Parallel::ThreadedActionCallback<ThreadedAction, ProxyType, Args...>>(
      cache, array_index, component_p, time, functions_to_check,
      std::forward<Args>(args)...);
}

namespace Actions {
/// \ingroup ComputationalDomainGroup
/// Check that functions of time are up-to-date.
///
/// Wait for all functions of time in `domain::Tags::FunctionsOfTime`
/// to be ready at `::Tags::Time`.  This ensures that the coordinates
/// can be safely accessed in later actions without first verifying
/// the state of the time-dependent maps.
template <size_t Dim>
struct CheckFunctionsOfTimeAreReady {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, ActionList /*meta*/,
      const ParallelComponent* component) {
    bool ready = false;

    if constexpr (Parallel::is_dg_element_collection_v<ParallelComponent>) {
      const auto element_location = static_cast<int>(
          Parallel::local_synchronous_action<
              Parallel::Actions::GetItemFromDistributedOject<
                  Parallel::Tags::ElementLocations<Dim>>>(
              Parallel::get_parallel_component<ParallelComponent>(cache))
              ->at(array_index));
      ASSERT(element_location == Parallel::my_node<int>(cache),
             "Expected to be running on node "
                 << Parallel::my_node<int>(cache)
                 << " but the record says it is on node " << element_location);
      const double time = db::get<::Tags::Time>(box);
      ready = domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime,
          Parallel::Actions::PerformAlgorithmOnElement<false>>(
          cache, element_location, component, time, std::nullopt, array_index,
          time);

      if ((not ready) and Parallel::get<control_system::Tags::Verbosity>(
                              cache) >= ::Verbosity::Debug) {
        Parallel::printf(
            "CheckFoTAction, %s: Functions of time are not ready at time "
            "%.16e. Registering threaded action callback\n",
            array_index, db::get<::Tags::Time>(box));
      }
    } else {
      ready = functions_of_time_are_ready_algorithm_callback<
          domain::Tags::FunctionsOfTime>(cache, array_index, component,
                                         db::get<::Tags::Time>(box));
    }

    return {ready ? Parallel::AlgorithmExecution::Continue
                  : Parallel::AlgorithmExecution::Retry,
            std::nullopt};
  }
};
}  // namespace Actions

/// \ingroup ComputationalDomainGroup
/// Dense-output postprocessor to check that functions of time are up-to-date.
///
/// Check that all functions of time in
/// `domain::Tags::FunctionsOfTime` are ready at `::Tags::Time`.  This
/// ensures that the coordinates can be safely accessed in later
/// actions without first verifying the state of the time-dependent
/// maps.  This postprocessor does not actually modify anything.
template <size_t Dim>
struct CheckFunctionsOfTimeAreReadyPostprocessor {
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<>;
  static void apply() {}

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ParallelComponent>
  static bool is_ready(
      const gsl::not_null<db::DataBox<DbTagsList>*> box,
      const gsl::not_null<tuples::TaggedTuple<InboxTags...>*> /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ParallelComponent* component) {
    if constexpr (Parallel::is_dg_element_collection_v<ParallelComponent>) {
      const auto element_location = static_cast<int>(
          Parallel::local_synchronous_action<
              Parallel::Actions::GetItemFromDistributedOject<
                  Parallel::Tags::ElementLocations<Dim>>>(
              Parallel::get_parallel_component<ParallelComponent>(cache))
              ->at(array_index));
      const double time = db::get<::Tags::Time>(*box);
      const bool ready =
          domain::functions_of_time_are_ready_threaded_action_callback<
              domain::Tags::FunctionsOfTime,
              Parallel::Actions::PerformAlgorithmOnElement<false>>(
              cache, element_location, component, time, std::nullopt,
              array_index, time);

      if ((not ready) and Parallel::get<control_system::Tags::Verbosity>(
                              cache) >= ::Verbosity::Debug) {
        Parallel::printf(
            "CheckFoTPostprocessor, %s: Functions of time are not ready at "
            "time "
            "%.16e. Registering threaded action callback\n",
            array_index, time);
      }

      return ready;
    } else {
      return functions_of_time_are_ready_algorithm_callback<
          domain::Tags::FunctionsOfTime>(cache, array_index, component,
                                         db::get<::Tags::Time>(*box));
    }
  }
};
}  // namespace domain
