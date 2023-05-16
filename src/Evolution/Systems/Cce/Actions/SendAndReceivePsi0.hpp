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
#include "Utilities/GetOutput.hpp"
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
      const TimeStepId& time_id = db::get<::Tags::TimeStepId>(box);
      if (Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(
              cache) >= ::Verbosity::Debug) {
        Parallel::printf("SendPsi0, t = %.16f: Sending Psi0.\n",
                         time_id.substep_time());
      }
      Parallel::receive_data<
          Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
          Parallel::get_parallel_component<
              typename Metavariables::gh_dg_element_array>(cache),
          time_id,
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
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    auto& inbox = tuples::get<
        Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
        inboxes);

    const Verbosity verbosity =
        Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(cache);

    const TimeStepId& time_id = db::get<::Tags::TimeStepId>(box);

    if (not Cce::is_outer_boundary_element<true>(make_not_null(&inbox), box,
                                                 cache)) {
      // if (verbosity >= ::Verbosity::Debug) {
      //   Parallel::printf(
      //       "ReceivePsi0 %s t = %.16f: Not a boundary element.
      //       Continuing.\n", get_output(array_index), time_id.substep_time());
      // }
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    if (inbox.empty()) {
      if (verbosity >= ::Verbosity::Debug) {
        std::stringstream ss{};
        bool first = true;
        tmpl::for_each<tmpl::list<InboxTags...>>(
            [&ss, &inboxes, &first](auto tag_v) {
              using tag = tmpl::type_from<decltype(tag_v)>;
              if (not first) {
                ss << ", ";
              } else {
                first = false;
              }
              ss << "(" << pretty_type::name<tag>() << ","
                 << tuples::get<tag>(inboxes).size() << ")";
            });

        Parallel::printf(
            "ReceivePsi0 %s t = %.16f: Inbox is empty. Waiting. All "
            "inboxes:\n %s\n",
            get_output(array_index), time_id.substep_time(), ss.str());
      }
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    if (verbosity >= ::Verbosity::Debug) {
      Parallel::printf("ReceivePsi0 %s t = %.16f: Inbox has time %.16f.\n",
                       get_output(array_index), time_id.substep_time(),
                       inbox.begin()->first.substep_time());
    }

    // Now we know we're on the boundary. Move the received data into the Inbox
    tmpl::for_each<typename Metavariables::ccm_psi0>(
        [&inbox, &box](auto tag_v) {
          using tag = typename decltype(tag_v)::type;
          db::mutate<tag>(
              make_not_null(&box),
              [&inbox](const gsl::not_null<typename tag::type*> destination) {
                *destination = get<tag>(inbox.begin()->second);
              });
        });

    inbox.clear();

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Cce::Actions
