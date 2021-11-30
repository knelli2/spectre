// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "ApparentHorizons/Strahlkorper.hpp"
#include "ApparentHorizons/Tags.hpp"
#include "ControlSystem/ApparentHorizons/Measurements.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Shape.hpp"
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
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Distorted;
}  // namespace Frame
/// \endcond

namespace control_system::Systems {
/*!
 * \brief Controls the \link domain::CoordinateMaps::TimeDependent::Shape Shape
 * \endlink map
 *
 * \details Controls the functions \f$ \lambda_{lm}(t) \f$ in the \link
 * domain::CoordinateMaps::TimeDependent::Shape Shape \endlink map to match the
 * shape of the excision sphere to the shape of the horizon.
 *
 * Requirements:
 * - This control system requires that there be at least one excision surface in
 *   the simulation
 * - Currently this control system can only be used with the \link
 *   control_system::ah::BothHorizons BothHorizons \endlink measurement
 * - Currently this control system can only be used with the \link
 *   control_system::ControlErrors::Shape Shape \endlink control error
 */
template <ah::HorizonLabel Horizon, size_t DerivOrder>
struct Shape : tt::ConformsTo<protocols::ControlSystem> {
  static constexpr size_t deriv_order = DerivOrder;

  static std::string name() {
    return Horizon == control_system::ah::HorizonLabel::AhA ? "ShapeA"
                                                            : "ShapeB";
  }

  using measurement = ah::BothHorizons;
  static_assert(tt::assert_conforms_to<measurement,
                                       control_system::protocols::Measurement>);

  using control_error = ControlErrors::Shape<Horizon>;
  static_assert(tt::assert_conforms_to<
                control_error, control_system::protocols::ControlError>);

  // tag goes in control component
  struct MeasurementQueue : db::SimpleTag {
    using type = LinkedMessageQueue<
        double, tmpl::list<QueueTags::Strahlkorper<Frame::Distorted>>>;
  };

  using simple_tags = tmpl::list<MeasurementQueue, Tags::ControlSystemName,
                                 Tags::ControlError<Shape>>;

  struct process_measurement {
    template <typename Submeasurement>
    using argument_tags =
        tmpl::list<StrahlkorperTags::Strahlkorper<Frame::Distorted>>;

    template <ah::HorizonLabel MeasureHorizon, typename Metavariables>
    static void apply(ah::BothHorizons::FindHorizon<MeasureHorizon> /*meta*/,
                      const Strahlkorper<Frame::Distorted>& strahlkorper,
                      Parallel::GlobalCache<Metavariables>& cache,
                      const LinkedMessageId<double>& measurement_id) {
      // The measurement event will call this for both horizons, but we only
      // need one of the horizons. So if it is called for the wrong horizon,
      // just do nothing.
      if constexpr (MeasureHorizon == Horizon) {
        auto& control_sys_proxy = Parallel::get_parallel_component<
            ControlComponent<Metavariables, Shape<Horizon, DerivOrder>>>(cache);

        Parallel::simple_action<::Actions::UpdateMessageQueue<
            QueueTags::Strahlkorper<Frame::Distorted>, MeasurementQueue,
            UpdateControlSystem<Shape>>>(control_sys_proxy, measurement_id,
                                         strahlkorper);
      } else {
        (void)strahlkorper;
        (void)cache;
        (void)measurement_id;
      }
    }
  };
};
}  // namespace control_system::Systems
