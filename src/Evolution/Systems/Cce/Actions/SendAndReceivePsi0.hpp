// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/Cce/IsBoundaryElement.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/CreateHasTypeAlias.hpp"

namespace Cce::Actions {
namespace detail {
CREATE_HAS_TYPE_ALIAS(gh_dg_element_array)
CREATE_HAS_TYPE_ALIAS_V(gh_dg_element_array)
}  // namespace detail

struct SendPsi0 {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    if constexpr (detail::has_gh_dg_element_array_v<Metavariables>) {
      Parallel::receive_data<
          Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
          Parallel::get_parallel_component<
              typename Metavariables::gh_dg_element_array>(cache),
          db::get<::Tags::TimeStepId>(box),
          db::get<::Tags::Variables<typename Metavariables::ccm_psi0>>(box),
          false);
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

template <typename Metavariables>
struct ReceivePsi0 {
  using inbox_tags = tmpl::list<
      Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>;
  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    auto& inbox = tuples::get<
        Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
        inboxes);

    if (not Cce::is_outer_boundary_element<true>(make_not_null(&inbox), box,
                                                 cache)) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    // Now we know we're on the boundary. Move the received data into the Inbox
    tmpl::for_each<typename Metavariables::ccm_psi0>(
        [&inbox, &box](auto tag_v) {
          using tag = typename decltype(tag_v)::type;
          db::mutate<tag>(
              make_not_null(&box),
              [&inbox](const gsl::not_null<typename tag::type*> destination,
                       const TimeStepId& time) {
                *destination = get<tag>(inbox[time]);
              },
              db::get<::Tags::TimeStepId>(box));
        });

    inbox.erase(db::get<::Tags::TimeStepId>(box));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Cce::Actions
