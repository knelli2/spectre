// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "ControlSystem/Component.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/TypeOfObservation.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Reduction.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/MakeString.hpp"

namespace control_system {
namespace detail {
struct WriterHelper {
  // Use funcl::AssertEqual for all because we aren't actually "reducing"
  // anything. These are just one time measurements.
  using ReductionData = Parallel::ReductionData<
      // Time
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
      // Lambda, dtLambda, d2tLambda
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
      // Q, dtQ, d2tQ
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
      Parallel::ReductionDatum<double, funcl::AssertEqual<>>>;

  const static inline std::vector<std::string> legend{
      "Time", "Lambda", "dtLambda", "d2tLambda", "Q", "dtQ", "d2tQ"};
};

template <typename ControlSystem>
struct Registration {
  template <typename ParallelComponent, typename DbTagsList,
            typename ArrayIndex>
  static std::pair<observers::TypeOfObservation, observers::ObservationKey>
  register_info(const db::DataBox<DbTagsList>& /*box*/,
                const ArrayIndex& /*array_index*/) {
    return {observers::TypeOfObservation::Reduction,
            observers::ObservationKey{ControlSystem::name()}};
  }
};
}  // namespace detail
template <typename ControlSystem, typename Metavariables>
void write_data_to_disk(
    const double time, Parallel::GlobalCache<Metavariables>& cache,
    const std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>&
        function_of_time,
    const std::array<DataVector, ControlSystem::deriv_order>& q_and_derivs,
    const DataVector& control_signal) {
  auto& observer_writer_proxy = Parallel::get_parallel_component<
      observers::ObserverWriter<Metavariables>>(cache);
  auto& control_component_proxy = Parallel::get_parallel_component<
      ControlComponent<Metavariables, ControlSystem>>(cache);

  constexpr size_t deriv_order = ControlSystem::deriv_order;
  std::array<DataVector, 3> function_at_current_time{};
  try {
    // If we are working with a QuaternionFunctionOfTime, we are actually
    // controlling a small angle so we want to write the angle to disk rather
    // than the quaternion. Since `omega_func_and_2_derivs` is not a virtual
    // function of the FunctionOfTime base class, we need to down cast the
    // original function
    const auto& quat_func_of_time = dynamic_cast<
        domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>&>(
        *function_of_time);
    function_at_current_time = quat_func_of_time.angle_func_and_2_derivs(time);
  } catch (const std::bad_cast& e) {
    // otherwise just call the usual `func_and_2_derivs` member.
    function_at_current_time = function_of_time->func_and_2_derivs(time);
  }

  // We have a different subfile for each component so loop over them
  const size_t num_components = function_at_current_time[0].size();
  for (size_t i = 0; i < num_components; ++i) {
    const std::string component_name = ControlSystem::component_name(i);
    // Currently all reduction data is written to the reduction file so preface
    // everything with ControlSystems/
    const std::string subfile_name{"/ControlSystems/" + ControlSystem::name() +
                                   component_name};
    const auto observation_id =
        observers::ObservationId(time, ControlSystem::name());
    const auto& legend = detail::WriterHelper::legend;

    Parallel::threaded_action<observers::ThreadedActions::WriteReductionData>(
        // Node 0 is always the writer
        observer_writer_proxy[0], observation_id,
        static_cast<size_t>(
            Parallel::my_node(*control_component_proxy.ckLocal())),
        subfile_name, legend,
        detail::WriterHelper::ReductionData{
            // clang-format off
            time,
            function_at_current_time[0][i],
            function_at_current_time[1][i],
            function_at_current_time[2][i],
            q_and_derivs[0][i],
            q_and_derivs[1][i],
            // d2tQ is -control_signal
            -control_signal[i]}  // clang-format on
    );
  }
}
}  // namespace control_system
