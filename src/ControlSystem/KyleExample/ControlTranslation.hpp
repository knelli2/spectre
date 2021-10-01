// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Translation.hpp"
#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "ControlSystem/KyleExample/TrackTranslation.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/UpdateControlSystem.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Actions/UpdateMessageQueue.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

#include <vector>
#include "Parallel/Printf.hpp"

namespace control_system {
template <size_t Dim>
struct Translation : tt::ConformsTo<protocols::ControlSystem> {
  static constexpr size_t deriv_order = Dim + 1;

  static std::string name() noexcept {
    return pretty_type::short_name<Translation<Dim>>();
  }
  using measurement = TrackTranslation;

  // tag goes in control component
  struct MeasurementQueue : db::SimpleTag {
    using type =
        LinkedMessageQueue<double, tmpl::list<QueueTags::MeasureTranslation>>;
  };

  using simple_tags =
      tmpl::list<MeasurementQueue, Tags::MeasureTranslationResult,
                 Tags::ControlSystemName>;

  struct process_measurement {
    template <typename Submeasurement>
    using argument_tags = tmpl::list<Tags::MeasureTranslationResult>;

    template <typename Metavariables>
    static void apply(SubTrackTranslation /*meta*/,
                      const std::vector<double>& measurement_result,
                      Parallel::GlobalCache<Metavariables>& cache,
                      const LinkedMessageId<double>& measurement_id) {
      auto& control_sys_proxy = Parallel::get_parallel_component<
          ControlComponent<Metavariables, Translation<Dim>>>(cache);
      // Parallel::printf("Inside process_measurement apply\n");
      Parallel::simple_action<::Actions::UpdateMessageQueue<
          // REPLACEFUNCTORHERE actually updates control system/funtions of time
          // print mesh velocity
          QueueTags::MeasureTranslation, MeasurementQueue,
          UpdateControlSystem<deriv_order, ControlErrors::Translation>>>(
          control_sys_proxy, measurement_id, measurement_result);
    }
  };
};
}  // namespace control_system
