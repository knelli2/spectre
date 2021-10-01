// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Parallel/Callback.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/StdHelpers.hpp"

#include "Parallel/Printf.hpp"

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
/// \endcond

namespace domain {
/// \ingroup ComputationalDomainGroup
/// Check that functions of time are up-to-date.
///
/// Check that functions of time in \p CacheTag with names in \p
/// functions_to_check are ready at time \p time.  If no names are
/// listed in \p functions_to_check, checks all functions in \p
/// CacheTag.  If any function is not ready, schedules a
/// `Parallel::PerformAlgorithmCallback` with the GlobalCache.
template <typename CacheTag, typename Metavariables, typename ArrayIndex,
          typename Component, size_t N = 0>
bool functions_of_time_are_ready(
    Parallel::GlobalCache<Metavariables>& cache, const ArrayIndex& array_index,
    const Component* /*meta*/, const double time,
    const std::array<std::string, N>& functions_to_check =
        std::array<std::string, 0>{}) noexcept {
  const auto& proxy =
      ::Parallel::get_parallel_component<Component>(cache)[array_index];

  return Parallel::mutable_cache_item_is_ready<CacheTag>(
      cache, [&functions_to_check, &proxy,
              &time](const std::unordered_map<
                     std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
                         functions_of_time) noexcept {
        ASSERT(alg::all_of(functions_to_check,
                           [&functions_of_time](
                               const std::string& function_to_check) noexcept {
                             return functions_of_time.count(
                                        function_to_check) == 1;
                           }),
               "Not all functions to check ("
                   << functions_to_check << ") are in the global cache ("
                   << keys_of(functions_of_time) << ")");
        for (const auto& [name, f_of_t] : functions_of_time) {
          if ((not functions_to_check.empty()) and
              alg::find(functions_to_check, name) == functions_to_check.end()) {
            continue;
          }
          const double expiration_time = f_of_t->time_bounds()[1];
          // Changed from >= to just > because a function of time is still valid
          // at its expiration time.
          if (time > expiration_time) {
            return std::unique_ptr<Parallel::Callback>(
                new Parallel::PerformAlgorithmCallback(proxy));
          }
        }
        return std::unique_ptr<Parallel::Callback>{};
      });
}

namespace Actions {
/// \ingroup ComputationalDomainGroup
/// Check that functions of time are up-to-date.
///
/// Wait for all functions of time in `domain::Tags::FunctionsOfTime`
/// to be ready at `::Tags::Time`.  This ensures that the coordinates
/// can be safely accessed in later actions without first verifying
/// the state of the time-dependent maps.
struct CheckFunctionsOfTimeAreReady {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTags>&&, Parallel::AlgorithmExecution> apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, ActionList /*meta*/,
      const ParallelComponent* component) noexcept {
    const bool ready =
        functions_of_time_are_ready<domain::Tags::FunctionsOfTime>(
            cache, array_index, component, db::get<::Tags::Time>(box));
    return {std::move(box), ready ? Parallel::AlgorithmExecution::Continue
                                  : Parallel::AlgorithmExecution::Retry};
  }
};
}  // namespace Actions
}  // namespace domain
