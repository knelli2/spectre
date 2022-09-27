// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "DataStructures/DataVector.hpp"
#include "Utilities/TMPL.hpp"

namespace control_system {
/*!
 * \ingroup ControlSystemGroup
 * \brief Construct the initial expiration times for functions of time that are
 * controlled by a control system
 *
 * The expiration times are constructed using inputs from control system
 * OptionHolders as an unordered map from the name of the function of time being
 * controlled to the expiration time. The expiration time is computed as
 * \f$\tau_\mathrm{exp} = \alpha_\mathrm{update} \tau_\mathrm{damp}\f$ where
 * \f$\alpha_\mathrm{update}\f$ is the update fraction supplied as input to the
 * Controller and \f$\tau_\mathrm{damp}\f$ is/are the damping timescales
 * supplied from the TimescaleTuner (\f$\tau_\mathrm{damp}\f$ is a DataVector
 * with as many components as the corresponding function of time, thus
 * \f$\tau_\mathrm{exp}\f$ will also be a DataVector of the same length).
 *
 * To protect against bad inputs, if the initial expiration time that is
 * calculated is smaller than the initial time step, then the expiration time is
 * simply set to the initial time step. However, the MeasurementTimescales have
 * the same protection so if this does happen, then something is most likely
 * wrong with your initial parameters for the control system.
 */
template <typename... OptionHolders>
std::unordered_map<std::string, double> initial_expiration_times(
    const double initial_time, const double initial_time_step,
    const OptionHolders&... option_holders) {
  std::unordered_map<std::string, double> initial_expiration_times{};

  [[maybe_unused]] const auto gather_initial_expiration_times =
      [&initial_time, &initial_time_step,
       &initial_expiration_times](const auto& option_holder) {
        const auto& controller = option_holder.controller;
        const std::string& name =
            std::decay_t<decltype(option_holder)>::control_system::name();
        const auto& tuner = option_holder.tuner;

        const double update_fraction = controller.get_update_fraction();
        const double curr_timescale = min(tuner.current_timescale());
        const double initial_expiration_time = update_fraction * curr_timescale;
        initial_expiration_times[name] =
            initial_time + std::max(initial_time_step, initial_expiration_time);
      };

  EXPAND_PACK_LEFT_TO_RIGHT(gather_initial_expiration_times(option_holders));

  return initial_expiration_times;
}

/*!
 * \ingroup ControlSystemGroup
 * \brief Calculate the next expiration time for the FunctionsOfTime.
 *
 * This is done by multiplying the new measurement timescale by the number of
 * measurements per update stored in the
 * control_system::Tags::MeasurementsPerUpdate tag and adding it to the current
 * time.
 */
double function_of_time_expiration_time(
    const double time, const DataVector& measurement_timescales,
    const int measurements_per_update);

/*!
 * \ingroup ControlSystemGroup
 * \brief Calculate the next expiration time for the MeasurementTimescales.
 *
 * First the function_of_time_expiration_time() is calculated. Then, half of the
 * minimum measurement timescale is subtracted off from that time and returned.
 * The reason for this is as follows:
 *
 * The last measurement for updating the control system happens at the function
 * of time expiration time \f$ t_{\mathrm{fot\ update}} \f$. Based on how dense
 * triggers are set up, which control_system::Trigger is a dense trigger, you
 * calculate the next trigger (measurement) time at the current measurement
 * time. However, at \f$ t_{\mathrm{fot\ update}} \f$ we don't know when the
 * next measurement time is because we haven't actually done the update and the
 * new measurement timescales are based on the updated damping timescales. So
 * since the measurement timescales are stored as 0th order
 * PiecewisePolynomials, if they expire at \f$ t_{\mathrm{fot\ update}} \f$ and
 * are also evaluated at \f$ t_{\mathrm{fot\ update}} \f$, they will give the
 * *old* measurement timescale and that will be used as the next measurement
 * after the update would happen.
 *
 * We don't want this. We want to wait until the control system has been updated
 * in order to choose the next measurement time. Thus, we make the measurement
 * timescales expire *earlier* than \f$ t_{\mathrm{fot\ update}} \f$ so that the
 * calculation of the next measurement timescale waits until the control system
 * (and the measurement timescales) have been updated. The fact that the
 * expiration time of the measurement timescales is half a measurement before
 * \f$ t_{\mathrm{fot\ update}} \f$ is just to ensure we are more than epsilon
 * before the last measurement and more than epsilon after the second to last
 * measurement (i.e. guaranteed to be between the last two measurements).
 */
double measurement_expiration_time(const double time,
                                   const DataVector& measurement_timescales,
                                   const int measurements_per_update);
}  // namespace control_system
