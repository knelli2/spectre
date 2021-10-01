// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "ControlSystem/Actions.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "ControlSystem/KyleExample/TrackTranslation.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Actions/UpdateMessageQueue.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

#include "Parallel/Printf.hpp"
#include <vector>

namespace control_system {
struct ControlTranslation : tt::ConformsTo<protocols::ControlSystem> {
  static std::string name() noexcept { return "Translation"; }
  using measurement = TrackTranslation;

  // tag goes in LinkedMessageQueue
  struct SubmeasureQueueTag {
    //using type = ::domain::Tags::MeshVelocity<2>::type;
    using type = std::vector<double>;
  };

  // tag goes in control component
  struct MeasurementQueue : db::SimpleTag {
    using type = LinkedMessageQueue<double, tmpl::list<SubmeasureQueueTag>>;
  };

  using simple_tags =
      tmpl::list<MeasurementQueue, ::Tags::MeasureTranslationResult,
                 Tags::ControlSystemName>;

  struct process_measurement {
    template <typename Submeasurement>
    using argument_tags = tmpl::list<::Tags::MeasureTranslationResult>;

    template <typename Metavariables>
    static void apply(SubTrackTranslation /*meta*/,
                      //const typename ::domain::Tags::MeshVelocity<2>::type&
                      const std::vector<double>&
                          measurement_result,
                      Parallel::GlobalCache<Metavariables>& cache,
                      const LinkedMessageId<double>& measurement_id) {
      auto& control_sys_proxy = Parallel::get_parallel_component<
          ControlComponent<Metavariables, ControlTranslation>>(cache);
      //Parallel::printf("Inside process_measurement apply\n");
      Parallel::simple_action<::Actions::UpdateMessageQueue<
          // REPLACEFUNCTORHERE actually updates control system/funtions of time
          // print mesh velocity
          SubmeasureQueueTag, MeasurementQueue, Actions::ControlMeshVelocity>>(
          control_sys_proxy, measurement_id, measurement_result);
    }
  };
};
}  // namespace control_system
