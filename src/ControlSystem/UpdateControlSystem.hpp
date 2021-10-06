// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace domain::Tags {
struct FunctionsOfTime;
}  // namespace domain::Tags

namespace control_system {
template <size_t DerivOrder, typename ControlError>
struct UpdateControlSystem {
  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename... TupleTags>
  static void apply(const gsl::not_null<db::DataBox<DbTags>*> box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time,
                    tuples::TaggedTuple<TupleTags...> data) {
    if constexpr (
        db::tag_is_retrievable_v<control_system::Tags::Averager<DerivOrder>,
                                 db::DataBox<DbTags>> and
        db::tag_is_retrievable_v<control_system::Tags::Controller<DerivOrder>,
                                 db::DataBox<DbTags>> and
        db::tag_is_retrievable_v<control_system::Tags::TimescaleTuner,
                                 db::DataBox<DbTags>> and
        db::tag_is_retrievable_v<::control_system::Tags::ControlSystemName,
                                 db::DataBox<DbTags>>) {
      const auto& functions_of_time =
          Parallel::get<::domain::Tags::FunctionsOfTime>(cache);
      const std::string& function_of_time_name =
          db::get<control_system::Tags::ControlSystemName>(*box);
      const auto& function_of_time =
          functions_of_time.at(function_of_time_name);

      // Compute control error from the template parameter
      const DataVector Q = ControlError::template apply<DerivOrder>(
          cache, time, function_of_time_name, data);

      // Get the averager, controller, and tuner from the box
      auto& averager =
          db::get_mutable_reference<control_system::Tags::Averager<DerivOrder>>(
              box);
      auto& controller = db::get_mutable_reference<
          control_system::Tags::Controller<DerivOrder>>(box);
      auto& tuner =
          db::get_mutable_reference<control_system::Tags::TimescaleTuner>(box);
      const DataVector& current_timescale = tuner.current_timescale();

      // Update the averager. We do this before we call controller.is_ready()
      // because we still want the averager to be up to date even if we aren't
      // updating at this time
      averager.update(time, Q, current_timescale);

      // Check if it is time to update
      if (not controller.is_ready(time)) {
        return;
      }

      // Get the averaged values of the control error and check if they are
      // valid
      const auto& opt_avg_values = averager(time);
      if (not opt_avg_values.has_value()) {
        return;
      }
      const double time_offset =
          averager.last_time_updated() - averager.average_time(time);
      const double time_offset_0th =
          averager.using_average_0th_deriv_of_q() ? time_offset : 0.0;

      // Calculate the control signal which will be used to update the highest
      // derivative of the FunctionOfTime
      const DataVector control_signal =
          controller(time, current_timescale, opt_avg_values.value(),
                     time_offset_0th, time_offset);

      const double current_expiration_time = function_of_time->time_bounds()[1];
      // Calculate the next expiration time based on the current one
      const double new_expiration_time =
          controller.next_expiration_time(current_expiration_time);

      // Actually update the FunctionOfTime
      Parallel::mutate<::domain::Tags::FunctionsOfTime, UpdateFunctionOfTime>(
          cache, function_of_time_name, current_expiration_time, control_signal,
          new_expiration_time);

      // Finally, update the timescales with the newly calculated control error
      // and its derivative
      std::array<DataVector, 2> updated_timescales{
          {opt_avg_values.value()[0], opt_avg_values.value()[1]}};
      tuner.update_timescale(updated_timescales);
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
