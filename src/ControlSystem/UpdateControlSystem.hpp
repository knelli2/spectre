// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"

#include "Parallel/Printf.hpp"
#include <ostream>

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
                    tuples::TaggedTuple<TupleTags...> data) noexcept {
    if constexpr (db::tag_is_retrievable_v<
                      control_system::Tags::Averager<DerivOrder>,
                      db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<
                      control_system::Tags::Controller<DerivOrder>,
                      db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<control_system::Tags::TimescaleTuner,
                                           db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<
                      ::control_system::Tags::ControlSystemName,
                      db::DataBox<DbTags>>) {
      // Parallel::printf("Hi from the good place\n");
      std::ostringstream os;
      os << "Time: " << time << "\n";

      const auto& functions_of_time =
          Parallel::get<::domain::Tags::FunctionsOfTime>(cache);
      const std::string& function_of_time_name =
          db::get<control_system::Tags::ControlSystemName>(*box);
      const auto& function_of_time =
          functions_of_time.at(function_of_time_name);

      const DataVector Q = ControlError::template apply<DerivOrder>(data);
      os << "Control error: " << Q << "\n";

      auto& averager = db::get_mutable_reference<
          control_system::Tags::Averager<DerivOrder>>(box);
      auto& controller =
          db::get<control_system::Tags::Controller<DerivOrder>>(*box);
      auto& tuner =
          db::get_mutable_reference<control_system::Tags::TimescaleTuner>(box);
      const DataVector& current_timescale = tuner.current_timescale();
      os << "  Current timescale: " << current_timescale << "\n";

      averager.update(time, Q, current_timescale);
      if (not controller.is_triggered(time)) {
        os << "  Controller not triggered.\n";
        if (db::get<control_system::Tags::PrintOutput>(*box)) {
          Parallel::printf(os.str());
        }
        return;
      }
      const auto& opt_avg_values = averager(time);
      os << "  Averager after update:\n" << averager.get_output() << "\n";
      if (not opt_avg_values.has_value()) {
        os << "  Not enough data to average\n";
        if (db::get<control_system::Tags::PrintOutput>(*box)) {
          Parallel::printf(os.str());
        }
        return;
      }
      // os << "  Averaged values: " << opt_avg_values.value() << "\n";
      const double time_offset =
          averager.last_time_updated() - averager.average_time(time);
      const double time_offset_0th =
          averager.using_average_0th_deriv_of_q() ? time_offset : 0.0;
      const DataVector control_signal =
          controller(time, current_timescale, opt_avg_values.value(),
                     time_offset_0th, time_offset);
      os << "  Control signal: " << control_signal << "\n";

      const double current_expiration_time = function_of_time->time_bounds()[1];
      // expiration time is next time I intend to update the fot
      const double new_expiration_time =
          controller.next_expiration_time(current_expiration_time);

      os << "  Updating FunctionOfTime with name '" << function_of_time_name
         << "' with:\n";
      os << "    curr_sim_time : " << time << "\n";
      os << "    curr_expr_time: " << current_expiration_time << "\n";
      os << "    control_signal: " << control_signal << "\n";
      os << "    new_expir_time: " << new_expiration_time << "\n";
      Parallel::mutate<::domain::Tags::FunctionsOfTime,
                       ::control_system::UpdateFunctionOfTime>(
          cache, function_of_time_name, current_expiration_time, control_signal,
          new_expiration_time);

      os << "  Updating timescales to be:\n";
      std::array<DataVector, 2> updated_timescales{
          {opt_avg_values.value()[0], opt_avg_values.value()[1]}};
      tuner.update_timescale(updated_timescales);
      os << "    " << tuner.current_timescale() << "\n";
      os << "  Next trigger will be at " << time + 0.3 * current_timescale[0]
         << "\n";
      if (db::get<control_system::Tags::PrintOutput>(*box)) {
        Parallel::printf(os.str());
      }
    } else {
      Parallel::printf("Hi from the bad place\n");
      const std::string& function_of_time_name =
          db::get<control_system::Tags::ControlSystemName>(*box);
      Parallel::printf(" Trying to find " + function_of_time_name + "\n");
      //(void)(box);
      (void)(time);
      (void)(data);
      // ERROR("Wrong DataBox.");
    }
    return;
  }
};
}  // namespace control_system
