// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <map>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "DataStructures/Variables.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Parallel/InboxInserters.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"

namespace Cce {
namespace ReceiveTags {

/// A receive tag for the data sent to the CCE evolution component from the CCE
/// boundary component
template <typename CommunicationTagList>
struct BoundaryData
    : Parallel::InboxInserters::Value<BoundaryData<CommunicationTagList>> {
  using temporal_id = TimeStepId;
  using type = std::unordered_map<temporal_id, Variables<CommunicationTagList>>;
};

/// A receive tag for the next time sent to the GH component from the CCE
/// evolution component
struct CcmNextTimeToGH : Parallel::InboxInserters::Value<CcmNextTimeToGH> {
  using temporal_id = TimeStepId;
  using type = std::map<temporal_id, TimeStepId>;
};

/// A receive tag for the element ID and next time sent to the CCM component
/// from all the GH boundary elements
///
/// Stored as a `std::unordered_map<ElementId, std::unordered_map<TimeStepId,
/// TimeStepId>>` where the first TimeStepId is the current time of the element
/// and the second `TimeStepId` is the next time of the element.
struct GhNextTimeToCcm {
  using temporal_id = TimeStepId;
  using mapped_type = std::pair<ElementId<3>, TimeStepId>;
  using type =
      std::unordered_map<ElementId<3>, std::pair<TimeStepId, TimeStepId>>;
  static void insert_into_inbox(const gsl::not_null<type*> inbox,
                                const temporal_id& received_temporal_id,
                                std::pair<ElementId<3>, TimeStepId>&& data) {
    (*inbox)[data.first] = std::make_pair(received_temporal_id, data.second);
  }
};
}  // namespace ReceiveTags
}  // namespace Cce
