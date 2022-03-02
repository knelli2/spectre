// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "ApparentHorizons/ObjectLabel.hpp"
#include "ApparentHorizons/Tags.hpp"
#include "ControlSystem/ApparentHorizons/Measurements.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Translation.hpp"
#include "ControlSystem/DataVectorHelpers.hpp"
#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Protocols/Measurement.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/UpdateControlSystem.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Actions/UpdateMessageQueue.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Grid;
}  // namespace Frame
/// \endcond

namespace control_system::Systems {
/*!
 * \brief Controls the 3D \link
 * domain::CoordinateMaps::TimeDependent::Translation Translation \endlink map
 *
 * \details Controls the function \f$ \vec{T}(t) \f$ in the \link
 * domain::CoordinateMaps::TimeDependent::Translation Translation \endlink map.
 *
 * Requirements:
 * - This control system requires that there be exactly two objects in the
 *   simulation
 * - Currently both these objects must be black holes
 * - Currently this control system can only be used with the \link
 *   control_system::ah::BothHorizons BothHorizons \endlink measurement
 * - Currently this control system can only be used with the \link
 *   control_system::ControlErrors::Translation Translation \endlink control
 *   error
 */
template <size_t DerivOrder>
struct Translation : tt::ConformsTo<protocols::ControlSystem> {
  static constexpr size_t deriv_order = DerivOrder;

  static std::string name() {
    return pretty_type::short_name<Translation<DerivOrder>>();
  }

  static std::string component_name(const size_t i) {
    return i == 0 ? "X" : (i == 1 ? "Y" : "Z");
  }

  using measurement = ah::BothHorizons;
  static_assert(tt::assert_conforms_to<measurement,
                                       control_system::protocols::Measurement>);

  using control_error = ControlErrors::Translation;
  static_assert(tt::assert_conforms_to<
                control_error, control_system::protocols::ControlError>);

  // tag goes in control component
  struct MeasurementQueue : db::SimpleTag {
    using type =
        LinkedMessageQueue<double,
                           tmpl::list<QueueTags::Center<::ah::ObjectLabel::A>,
                                      QueueTags::Center<::ah::ObjectLabel::B>>>;
  };

  using simple_tags = tmpl::list<MeasurementQueue, Tags::ControlSystemName,
                                 Tags::ControlError<Translation>>;

  struct process_measurement {
    template <typename Submeasurement>
    using argument_tags =
        tmpl::list<StrahlkorperTags::Strahlkorper<Frame::Grid>>;

    template <::ah::ObjectLabel Horizon, typename Metavariables>
    static void apply(ah::BothHorizons::FindHorizon<Horizon> /*meta*/,
                      const Strahlkorper<Frame::Grid>& strahlkorper,
                      Parallel::GlobalCache<Metavariables>& cache,
                      const LinkedMessageId<double>& measurement_id) {
      auto& control_sys_proxy = Parallel::get_parallel_component<
          ControlComponent<Metavariables, Translation<DerivOrder>>>(cache);

      const DataVector center =
          array_to_datavector(strahlkorper.physical_center());

      Parallel::simple_action<::Actions::UpdateMessageQueue<
          QueueTags::Center<Horizon>, MeasurementQueue,
          UpdateControlSystem<Translation>>>(control_sys_proxy, measurement_id,
                                             center);
    }
  };
};
}  // namespace control_system::Systems
