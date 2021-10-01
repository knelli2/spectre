// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <ostream>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Printf.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Tags {
struct Time;
}  // namespace Tags

namespace domain::Tags {
struct FunctionsOfTime;
}  // namespace domain::Tags

namespace control_system {
namespace Actions {
struct ControlMeshVelocity {
  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename SubmeasurementTag>
  static void apply(const gsl::not_null<db::DataBox<DbTags>*> box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time,
                    tuples::TaggedTuple<SubmeasurementTag> data) noexcept {
    constexpr size_t Dim = 1;
    if constexpr (db::tag_is_retrievable_v<::Tags::Averager<Dim + 1>,
                                           db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<::Tags::Controller<Dim + 1>,
                                           db::DataBox<DbTags>> and
                  db::tag_is_retrievable_v<::Tags::TimescaleTuner,
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

      const auto& integral_vector = get<SubmeasurementTag>(data);
      DataVector target_value{Dim, 0.0};
      // integral(psi^2 * x_i) / integral(psi^2)
      // basically a CoM
      for (size_t i = 1; i < Dim + 1; i++) {
        target_value[i - 1] = integral_vector[i] / integral_vector[0];
      }

      const double& center_of_interval = integral_vector[Dim + 1];
      os << "Center of interval: " << center_of_interval << "\n";
      const DataVector Q = target_value - center_of_interval;
      os << "Control error: " << Q << "\n";

      auto& averager =
          db::get_mutable_reference<::Tags::Averager<Dim + 1>>(box);
      auto& controller = db::get<::Tags::Controller<Dim + 1>>(*box);
      auto& tuner = db::get_mutable_reference<::Tags::TimescaleTuner>(box);
      const DataVector& current_timescale = tuner.current_timescale();
      os << "  Current timescale: " << current_timescale << "\n";

      averager.update(time, Q, current_timescale);
      if (not controller.is_triggered(time)) {
        os << "  Controller not triggered.\n";
        if (db::get<::Tags::PrintOutput>(*box)) {
          Parallel::printf(os.str());
        }
        return;
      }
      const auto& opt_avg_values = averager(time);
      os << "  Averager after update:\n" << averager.get_output() << "\n";
      if (not opt_avg_values.has_value()) {
        os << "  Not enough data to average\n";
        if (db::get<::Tags::PrintOutput>(*box)) {
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
      if (db::get<::Tags::PrintOutput>(*box)) {
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

template <typename ControlSystems>
struct ReduceToRunCallbacks {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(const db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const LinkedMessageId<double>& measurement_id,
                    const std::vector<double>& integrals,
                    const std::vector<DataVector>& vec_coords) noexcept {
    // const typename ::Tags::MeasureTranslationResult::type&
    //    data_from_element) noexcept {
    std::ostringstream os;
    // os << "Reduction finished! Inside ReduceToRunCallbacks.\n";
    // Parallel::printf(os.str());
    double coord_min = std::numeric_limits<double>::max();
    double coord_max = std::numeric_limits<double>::min();
    for (auto& dv : vec_coords) {
      if (min(dv) < coord_min) {
        coord_min = min(dv);
      }
      if (max(dv) > coord_max) {
        coord_max = max(dv);
      }
    }
    const double center_coord = 0.5 * (coord_max + coord_min);
    std::vector<double> integrals2 = integrals;
    integrals2.push_back(center_coord);
    const auto temp_box =
        db::create<db::AddSimpleTags<::Tags::MeasureTranslationResult>>(
            integrals2);
    //        data_from_element);
    control_system::RunCallbacks<control_system::SubTrackTranslation,
                                 ControlSystems>::apply(temp_box, cache,
                                                        measurement_id);
  }
};
}  // namespace Actions
}  // namespace control_system
