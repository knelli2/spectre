// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/Variables.hpp"
#include "Parallel/InboxInserters.hpp"
#include "Time/TimeStepId.hpp"

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
struct CcmNextTimeToGh : Parallel::InboxInserters::Value<CcmNextTimeToGh> {
  using temporal_id = TimeStepId;
  using type = std::map<temporal_id, double>;
};

}  // namespace ReceiveTags
}  // namespace Cce
