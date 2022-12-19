// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace Cce {
namespace Actions {

template <typename Metavariables>
struct CcmNextTime;

template <typename CceComponent>
struct SendNextTimeToGh {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    double next_time =
        db::get<::Tags::Next<::Tags::TimeStepId>>(box).substep_time().value();

    if constexpr (tmpl::size<tmpl::filter<
                      typename Metavariables::component_list,
                      tt::is_a<::DgElementArray, tmpl::_1>>>::value == 1) {
      Parallel::receive_data<Cce::ReceiveTags::CcmNextTimeToGh>(
          Parallel::get_parallel_component<
              typename Metavariables::gh_dg_element_array>(cache),
          db::get<::Tags::TimeStepId>(box), next_time, false);
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

} // namespace Actions
} // namespace Cce
