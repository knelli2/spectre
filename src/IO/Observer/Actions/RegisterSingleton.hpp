// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/Observer/Actions/ObserverRegistration.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/TypeOfObservation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/TaggedTuple.hpp"

#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/PrettyType.hpp"

namespace observers::Actions {
/*!
 * \brief Registers a singleton with the ObserverWriter.
 *
 * The singleton that observes is expected to call WriteReductionData
 * on node 0.
 */
template <typename RegisterHelper>
struct RegisterSingletonWithObserverWriter {
  template <typename DbTagList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagList>&&, bool> apply(
      db::DataBox<DbTagList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto [type_of_observation, observation_key] =
        RegisterHelper::template register_info<ParallelComponent>(box,
                                                                  array_index);

    switch (type_of_observation) {
      case TypeOfObservation::Reduction:
        break;
      case TypeOfObservation::Volume:
        ERROR(
            "Registering volume observations is not supported for singletons. "
            "The TypeOfObservation should be 'Reduction'.");
      default:
        ERROR(
            "Registering an unknown TypeOfObservation. It should be "
            "'Reduction' for singleton.");
    };

    // We call only on node 0; the observation call will occur only
    // on node 0.
    auto& my_proxy = Parallel::get_parallel_component<ParallelComponent>(cache);
    auto* my_ptr = my_proxy.ckLocal();
    ASSERT(my_ptr != nullptr,
           "Inside RegisterSingleton. The %s.ckLocal() gave a null pointer and "
           "it shouldn't have! What "
           "have you done! Oh the humanity..."
               << pretty_type::get_name<ParallelComponent>());
    Parallel::simple_action<Actions::RegisterReductionNodeWithWritingNode>(
        Parallel::get_parallel_component<
            observers::ObserverWriter<Metavariables>>(cache)[0],
        observation_key,
        static_cast<size_t>(Parallel::my_node(*my_proxy.ckLocal())));
    return {std::move(box), true};
  }
};
}  // namespace observers::Actions
