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
#include "ControlSystem/ControlErrors/Rotation.hpp"
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
 * \brief Controls the 3D \link domain::CoordinateMaps::TimeDependent::Rotation
 * Rotation \endlink map
 *
 * \details Indirectly controls the quaternion in the 3D Rotation coordinate
 * map. It does this by directly controlling the internal PiecewisePolynomial in
 * a \link domain::FunctionsOfTime::QuaternionFunctionOfTime
 * QuaternionFunctionOfTime \endlink, which represents the angle the system has
 * rotated through about each axis.
 *
 * Requirements:
 * - The function of time this control system controls must be a
 *   QuaternionFunctionOfTime.
 * - This control system requires that there be exactly two objects in the
 *   simulation
 * - Currently both these objects must be black holes
 * - Currently this control system can only be used with the \link
 *   control_system::ah::BothHorizons BothHorizons \endlink measurement
 * - Currently this control system can only be used with the \link
 *   control_system::ControlErrors::Rotation Rotation \endlink control error
 */
template <size_t DerivOrder>
struct Rotation : tt::ConformsTo<protocols::ControlSystem> {
  static constexpr size_t deriv_order = DerivOrder;

  static std::string name() {
    return pretty_type::short_name<Rotation<DerivOrder>>();
  }

  static std::string component_name(const size_t i) {
    return i == 0 ? "X" : (i == 1 ? "Y" : "Z");
  }

  using measurement = ah::BothHorizons;
  static_assert(tt::assert_conforms_to<measurement,
                                       control_system::protocols::Measurement>);

  using control_error = ControlErrors::Rotation;
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
                                 Tags::ControlError<Rotation>>;

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
          ControlComponent<Metavariables, Rotation<DerivOrder>>>(cache);

      const DataVector center =
          array_to_datavector(strahlkorper.physical_center());

      Parallel::simple_action<::Actions::UpdateMessageQueue<
          QueueTags::Center<Horizon>, MeasurementQueue,
          UpdateControlSystem<Rotation>>>(control_sys_proxy, measurement_id,
                                          center);
    }
  };
};
}  // namespace control_system::Systems
