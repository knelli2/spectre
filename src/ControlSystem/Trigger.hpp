// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "ControlSystem/Metafunctions.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/FunctionsOfTimeAreReady.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace Tags {
struct Time;
}  // namespace Tags
namespace control_system::Tags {
struct MeasurementTimescales;
}  // namespace control_system::Tags
/// \endcond

namespace control_system {
/// \ingroup ControlSystemsGroup
/// \ingroup EventsAndTriggersGroup
/// Trigger for control system measurements.
///
/// This trigger is only intended to be used with the
/// `control_system::Event` event.  A specialization of this trigger
/// will be created during control system initialization for each
/// unique \ref control_system::protocols::Measurement "measurement".
///
/// These triggers must be added to the \ref
/// Options::protocols::FactoryCreation "factory_creation" struct in
/// the metavariables, even though they cannot be created from the
/// input file.  The `control_system::control_system_triggers`
/// metafunction provides the list of triggers to include.
template <typename ControlSystems>
class Trigger : public DenseTrigger {
  static_assert(tmpl::size<ControlSystems>::value > 0);
  using measurement = typename tmpl::front<ControlSystems>::measurement;
  static_assert(tmpl::all<ControlSystems,
                          std::is_same<metafunctions::measurement<tmpl::_1>,
                                       tmpl::pin<measurement>>>::value);

 public:
  /// \cond
  // LCOV_EXCL_START
  explicit Trigger(CkMigrateMessage* const msg) : DenseTrigger(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Trigger);  // NOLINT
  // LCOV_EXCL_STOP
  /// \endcond

  // This trigger is created during control system initialization, not
  // from the input file.
  static constexpr bool factory_creatable = false;
  Trigger() = default;

  using is_triggered_argument_tags =
      tmpl::list<::Tags::Time, control_system::Tags::MeasurementTimescales,
                 domain::Tags::FunctionsOfTime>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  std::optional<bool> is_triggered(
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const Component* /*component*/,
      const double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          measurement_timescales,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) {
    if (UNLIKELY(not next_trigger_.has_value())) {
      // First call

      // This will happen if an executable has control systems, but
      // all functions of time were overriden by ones read in from a
      // file. So there is no need to trigger control systems.  Since
      // we only enter this branch on the first call to the trigger,
      // this is the initial time so we can assume the
      // measurement_timescales are ready.
      if (next_measurement(time, measurement_timescales, functions_of_time) ==
          std::numeric_limits<double>::infinity()) {
        next_trigger_ = std::numeric_limits<double>::infinity();
      } else {
        next_trigger_ = time;
      }
    }

    return time == *next_trigger_;
  }

  using next_check_time_argument_tags =
      tmpl::list<::Tags::Time, control_system::Tags::MeasurementTimescales,
                 domain::Tags::FunctionsOfTime>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  std::optional<double> next_check_time(
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const Component* component,
      const double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          measurement_timescales,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) {
    // At least one control system is active
    const bool is_ready = tmpl::as_pack<ControlSystems>(
        [&array_index, &cache, &component, &time](auto... control_systems) {
          return domain::functions_of_time_are_ready<
              control_system::Tags::MeasurementTimescales>(
              cache, array_index, component, time,
              std::array{
                  tmpl::type_from<decltype(control_systems)>::name()...});
        });
    if (not is_ready) {
      return std::nullopt;
    }

    const bool triggered = time == *next_trigger_;
    if (triggered) {
      *next_trigger_ =
          next_measurement(time, measurement_timescales, functions_of_time);
    }
    return *next_trigger_;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    DenseTrigger::pup(p);
    p | next_trigger_;
  }

 private:
  static double next_measurement(
      const double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          measurement_timescales,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) {
    const double min_measure_time = tmpl::as_pack<ControlSystems>(
        [&measurement_timescales, &time](auto... control_systems) {
          return std::min(
              {min(measurement_timescales
                       .at(tmpl::type_from<decltype(control_systems)>::name())
                       ->func(time)[0])...});
        });
    if (min_measure_time == std::numeric_limits<double>::infinity()) {
      return min_measure_time;
    }
    const double calculated_next_measurement = time + min_measure_time;

    double min_fot_expr_time = std::numeric_limits<double>::infinity();
    for (const auto& [name, fot] : functions_of_time) {
      (void)name;
      min_fot_expr_time = std::min(min_fot_expr_time, fot->time_bounds()[1]);
    }

    // If the expiration time is infinity, don't bother checking if it's equal
    // to the min measure time. We already returned above if it is infinity
    if (min_fot_expr_time == std::numeric_limits<double>::infinity()) {
      return calculated_next_measurement;
    }

    double next_measurement_time;
    // Just take the average of the two things we are comparing. They'll always
    // be super close
    const double scale =
        0.5 * (calculated_next_measurement + min_fot_expr_time);
    // To avoid roundoff differences between expiration times and when
    // measurements are supposed to happen, if the next measurement is within
    // roundoff of the current expiration time, use the exact value of the
    // expiration time as the measurement time.
    if (equal_within_roundoff(calculated_next_measurement, min_fot_expr_time,
                              std::numeric_limits<double>::epsilon() * 100.0,
                              scale)) {
      next_measurement_time = min_fot_expr_time;
    } else {
      next_measurement_time = calculated_next_measurement;
    }

    return next_measurement_time;
  }

  std::optional<double> next_trigger_{};
};

/// \cond
template <typename ControlSystems>
PUP::able::PUP_ID Trigger<ControlSystems>::my_PUP_ID = 0;  // NOLINT
/// \endcond

// This metafunction is tested in Test_EventTriggerMetafunctions.cpp

/// \ingroup ControlSystemGroup
/// The list of triggers needed for measurements for a list of control
/// systems.
template <typename ControlSystems>
using control_system_triggers = tmpl::transform<
    metafunctions::measurements_t<ControlSystems>,
    tmpl::bind<Trigger, metafunctions::control_systems_with_measurement<
                            tmpl::pin<ControlSystems>, tmpl::_1>>>;
}  // namespace control_system
