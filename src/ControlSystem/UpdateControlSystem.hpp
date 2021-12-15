// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"

#include <ostream>
#include "Parallel/Printf.hpp"

namespace domain::Tags {
struct FunctionsOfTime;
}  // namespace domain::Tags

namespace control_system {
/*!
 * \brief Functor for updating control systems when they are ready.
 *
 * \details The apply operator of this struct is meant to be called by the
 * UpdateMessageQueue action once an entire measurement has been made.
 *
 * Requires a few tags to be in the DataBox of the ControlComponent that this
 * is running on:
 * - \link Tags::Averager Averager \endlink
 * - \link Tags::Controller Controller \endlink
 * - \link Tags::TimescaleTuner TimescaleTuner \endlink
 * - \link Tags::ControlSystemName ControlSystemName \endlink
 * - \link Tags::ControlError ControlError \endlink
 *
 * If these are not in the DataBox, an error will occur.
 *
 * The algorithm to determine whether or not to update the functions of time is
 * as follows:
 * 1. Determine if we need to update now. This is done by checking if the next
 *    measurement scheduled to happen will be after the current expiration time
 *    of the function of time we are controlling. If it is we need to update
 *    now, otherwise the functions of time will be expired the next time we
 *    measure.
 * 2. Determine if we should store the current measurement. Not all control
 *    systems are on the same timescale; some are on shorter timescales (like
 *    characteristic speed control), some are on longer timescales (like
 *    Rotation). The ones on longer timescales don't need to record measurements
 *    as often as the ones on shorter timescales. If they did, this could
 *    introduce unnecessary noise into the measurements. It also causes the
 *    control errors to be calculated based off only measurements from the end
 *    of the update interval. This is not what we want. We want measurements
 *    throughout the udpate interval. If we don't need to store the current
 *    measurement and we don't need to update now, end here.
 * 3. Calculate the control error and update the averager (store the current
 *    measurement).
 * 4. Determine if we need to update. We only want to update at the very end of
 *    the update interval, not sooner. If we don't need to update, end here.
 * 5. Compute control signal using the control error and its derivatives.
 * 6. Determine the new expiration time.
 * 7. Update the function of time.
 * 8. Update the damping timescale using the control error and one derivative.
 * 9. Determine the new measurement timescale for this function of time using
 *    the new damping timescale.
 * 10. Update the measurement timescale.
 */
template <typename ControlSystem>
struct UpdateControlSystem {
  static constexpr size_t deriv_order = ControlSystem::deriv_order;

  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename... TupleTags>
  static void apply(const gsl::not_null<db::DataBox<DbTags>*> box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time,
                    tuples::TaggedTuple<TupleTags...> data) {
    if constexpr (db::tag_is_retrievable_v<
                      control_system::Tags::Averager<deriv_order - 1>,
                      db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<
                      control_system::Tags::Controller<deriv_order>,
                      db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<control_system::Tags::TimescaleTuner,
                                           db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<
                      control_system::Tags::ControlError<ControlSystem>,
                      db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<
                      ::control_system::Tags::ControlSystemName,
                      db::DataBox<DbTags>>) {
      // Begin step 1
      const auto& functions_of_time =
          Parallel::get<::domain::Tags::FunctionsOfTime>(cache);
      const auto& measurement_timescales =
          Parallel::get<Tags::MeasurementTimescales>(cache);
      const std::string& function_of_time_name =
          db::get<control_system::Tags::ControlSystemName>(*box);
      const auto& function_of_time =
          functions_of_time.at(function_of_time_name);

      const double current_expiration_time = function_of_time->time_bounds()[1];
      // We need the min over all functions of time because that's when the next
      // measurement will be
      double current_min_time_between_measurements{
          std::numeric_limits<double>::max()};
      for (const auto& [name, measurement_timescale_fot] :
           measurement_timescales) {
        current_min_time_between_measurements =
            std::min(current_min_time_between_measurements,
                     min(measurement_timescale_fot->func(time)[0]));
      }

      // If the next time we are going to measure is after the current
      // expiration time, we have to update things now
      const bool need_to_update_now =
          current_expiration_time <=
          time + current_min_time_between_measurements;
      // End step 1

      // Begin step 2
      // Get the averager, controller, tuner, and control error from the box
      auto& averager = db::get_mutable_reference<
          control_system::Tags::Averager<deriv_order - 1>>(box);
      auto& controller = db::get_mutable_reference<
          control_system::Tags::Controller<deriv_order>>(box);
      auto& tuner =
          db::get_mutable_reference<control_system::Tags::TimescaleTuner>(box);
      auto& control_error = db::get_mutable_reference<
          control_system::Tags::ControlError<ControlSystem>>(box);
      const DataVector& current_timescale = tuner.current_timescale();

      // Check if we actually need to use the measurement. If a control system
      // is on a longer timescale, we don't need to store frequent measurements
      if (not(averager.is_ready(time) or need_to_update_now)) {
        return;
      }
      // End step 2

      // Begin step 3
      // Compute control error
      const DataVector Q =
          control_error(cache, time, function_of_time_name, data);

      // Update the averager. We do this before we call controller.is_ready()
      // because we still want the averager to be up to date even if we aren't
      // updating at this time
      averager.update(time, Q, current_timescale);
      // End step 3

      // Begin step 4
      // Check if it is time to update
      if (not(controller.is_ready(time) or need_to_update_now)) {
        return;
      }
      // End step 4

      // Begin step 5
      // Get the averaged values of the control error and its derivatives
      const auto& opt_avg_values = averager(time);

      const double time_offset =
          averager.last_time_updated() - averager.average_time(time);
      const double time_offset_0th =
          averager.using_average_0th_deriv_of_q() ? time_offset : 0.0;

      // Calculate the control signal which will be used to update the highest
      // derivative of the FunctionOfTime
      const DataVector control_signal =
          controller(time, current_timescale, opt_avg_values.value(),
                     time_offset_0th, time_offset);
      // End step 5

      // Begin step 6
      // Calculate the next expiration time based on the current one
      const double new_expiration_time =
          controller.next_expiration_time(current_expiration_time);
      // End step 6

      // Begin step 7
      // Actually update the FunctionOfTime
      Parallel::mutate<::domain::Tags::FunctionsOfTime, UpdateFunctionOfTime>(
          cache, function_of_time_name, current_expiration_time, control_signal,
          new_expiration_time);
      // End step 7

      // Begin step 8
      // Update the damping timescales with the newly calculated control error
      // and its derivative
      std::array<DataVector, 2> q_and_dtq{
          {(*opt_avg_values)[0], {(*opt_avg_values)[0].size(), 0.0}}};
      if constexpr (deriv_order > 1) {
        q_and_dtq[1] = (*opt_avg_values)[1];
      }
      tuner.update_timescale(q_and_dtq);
      // End step 8

      // Begin step 9
      // Calculate new measurement timescales with updated damping timescales
      const DataVector new_measurement_timescale =
          Tags::detail::calculate_measurement_timescales(controller, tuner);
      // End step 9

      // Begin step 10
      // Update the measurement timescales
      Parallel::mutate<Tags::MeasurementTimescales, UpdateFunctionOfTime>(
          cache, function_of_time_name, current_expiration_time,
          new_measurement_timescale, new_expiration_time);

      // Now that the measurement timescales have been updated, tell the
      // averager when to expect the next measurement
      averager.assign_time_between_measurements(min(new_measurement_timescale));
      // End step 10
    } else {
      (void)(box);
      (void)(cache);
      (void)(time);
      (void)(data);
      ERROR(
          "Wrong DataBox. Is this the DataBox for the correct control "
          "component?");
    }
  }
};
}  // namespace control_system
