// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <memory>
#include <string>
#include <unordered_map>

#include "ControlSystem/Controller.hpp"
#include "ControlSystem/InitialExpirationTimes.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Time/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
template <class Metavariables, typename ControlSystem>
struct ControlComponent;
/// \endcond

namespace control_system::Tags {
namespace detail {
template <size_t DerivOrder>
DataVector calculate_measurement_timescales(
    const ::Controller<DerivOrder>& controller, const ::TimescaleTuner& tuner) {
  const DataVector timescale = tuner.current_timescale();
  const double update_fraction = controller.get_update_fraction();
  const double frac = 1.0 / static_cast<double>(DerivOrder + 1);
  const DataVector measure_timescale_1 = timescale * update_fraction * frac;
  const DataVector measure_timescale_2 =
      tuner.current_timescale() * controller.get_update_fraction() *
      (1.0 / static_cast<double>(DerivOrder + 1));
  Parallel::printf(
      "Inside calculate_measurement_timescales:\n"
      " Deriv order: %d\n"
      " update fraction: %.17f\n"
      " damp timescale: %s\n"
      " frac: %.17f\n"
      " measure timescale 1: %s\n"
      " measure timescale 2: %s\n",
      DerivOrder, update_fraction, timescale, frac, measure_timescale_1,
      measure_timescale_2);
  return tuner.current_timescale() * controller.get_update_fraction() *
         (1.0 / static_cast<double>(DerivOrder + 1));
}
}  // namespace detail

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// \brief The measurement timescales associated with
/// domain::Tags::FunctionsOfTime.
///
/// Each function of time associated with a control system has a corresponding
/// set of timescales here, represented as `PiecewisePolynomial<0>` with the
/// same components as the function itself.
struct MeasurementTimescales : db::SimpleTag {
  using type = std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;

  static constexpr bool pass_metavariables = true;

  template <typename Metavariables>
  using option_tags = tmpl::flatten<
      tmpl::list<::OptionTags::InitialTime, ::OptionTags::InitialTimeStep,
                 control_system::inputs<tmpl::transform<
                     tmpl::filter<typename Metavariables::component_list,
                                  tt::is_a_lambda<ControlComponent, tmpl::_1>>,
                     tmpl::bind<tmpl::back, tmpl::_1>>>>>;

  template <typename Metavariables, typename... OptionHolders>
  static type create_from_options(const double initial_time,
                                  const double initial_time_step,
                                  const OptionHolders&... option_holders) {
    std::unordered_map<std::string,
                       std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
        timescales{};
    const auto initial_expiration_times =
        control_system::initial_expiration_times(
            initial_time, initial_time_step, option_holders...);

    [[maybe_unused]] const auto calculate_measurement_timescales =
        [&timescales, &initial_time, &initial_time_step,
         &initial_expiration_times](const auto& option_holder) {
          // This check is intentionally inside the lambda so that it will not
          // trigger for domains without control systems.
          if (initial_time_step <= 0.0) {
            ERROR(
                "Control systems can only be used in forward-in-time "
                "evolutions.");
          }

          const auto& controller = option_holder.controller;
          const std::string& name =
              std::decay_t<decltype(option_holder)>::control_system::name();
          const auto& tuner = option_holder.tuner;

          DataVector measurement_timescales =
              detail::calculate_measurement_timescales(controller, tuner);
          // At a minimum, we can only measure once a time step with GTS.
          for (size_t i = 0; i < measurement_timescales.size(); i++) {
            measurement_timescales[i] =
                std::max(initial_time_step, measurement_timescales[i]);
          }

          timescales.emplace(
              name,
              std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
                  initial_time, std::array{std::move(measurement_timescales)},
                  initial_expiration_times.at(name)));
        };

    EXPAND_PACK_LEFT_TO_RIGHT(calculate_measurement_timescales(option_holders));

    return timescales;
  }
};
}  // namespace control_system::Tags
