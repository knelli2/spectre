// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/MemoryMonitor/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace Initialization {
namespace Actions {
/*!
 * \brief Initialize the MemoryMonitor parallel component
 *
 * \details Adds a mem_monitor::Tags::MemoryHolder tag to the DataBox for every
 * parallel component in the `component_list` type alias in the metavariables.
 * Also adds a tag for the GlobalCache and MutableGlobalCache because we may
 * want to keep track of their size as well. All tags are constructed by default
 * even if we aren't monitoring every component.
 */
template <typename Metavariables>
struct InitializeMemoryMonitor {
  using component_list = typename Metavariables::component_list;

  using simple_tags = tmpl::push_back<
      tmpl::transform<component_list,
                      tmpl::bind<mem_monitor::Tags::MemoryHolder, tmpl::_1>>,
      mem_monitor::Tags::MemoryHolder<Parallel::GlobalCache<Metavariables>>,
      mem_monitor::Tags::MemoryHolder<
          Parallel::MutableGlobalCache<Metavariables>>>;

  using compute_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/, ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace Actions
}  // namespace Initialization
