// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"

namespace Cce {
namespace Actions {

template <typename Metavariables>
struct ReceiveCcmNextTime {
  using inbox_tags = tmpl::list<Cce::ReceiveTags::CcmNextTimeToGh>;
  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    auto& inbox = tuples::get<Cce::ReceiveTags::CcmNextTimeToGh>(inboxes);

    const auto& element = db::get<domain::Tags::Element<3_st>>(box);
    if (element.external_boundaries().size() == 0_st) {
      inbox.clear();
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    double next_time = std::numeric_limits<double>::signaling_NaN();
    if (inbox.size() > 0) {
      next_time = inbox.begin()->second;
    } else {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    db::mutate<Cce::Tags::NextCcmTime>(
        make_not_null(&box),
        [&next_time](const gsl::not_null<double*> next_ccm_time) {
            *next_ccm_time = next_time;
        });

    inbox.erase(inbox.begin());
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace Actions
}  // namespace Cce
