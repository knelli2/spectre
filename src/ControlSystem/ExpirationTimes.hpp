// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "ControlSystem/CalculateMeasurementTimescales.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Utilities/TMPL.hpp"

namespace control_system {
/*!
 * \ingroup ControlSystemGroup
 * \brief Calculate the next expiration time for the FunctionsOfTime.
 *
 * This is done by adding to the current time one old minimum measurement
 * timescale and measurements per update times the new minimum measurement
 * timescale. It is done this way because we update the functions of time one
 * (old) measurement before they actually expire.
 *
 * The choice of having the functions of time expire exactly one old measurement
 * after they are updated is arbitrary. They could expire any time between the
 * update time and one old measurement after the update. This decision was made
 * to minimize time spent waiting for the functions of time to be valid.
 *
 * Since functions of time are valid at their expiration time, we are actually
 * able to do the next measurement if the expiration time is at that
 * measurement. Thus we delay any potential waiting that may happen until the
 * subsequent measurement (and by that time, most, if not all, functions of time
 * should have been updated because an entire horizon find happened in the
 * meantime). If the expiration time was earlier than the next measurement, we'd
 * have to pause the Algorithm on the elements and wait until all the functions
 * of time have been updated.
 */
double function_of_time_expiration_time(
    const double time, const DataVector& old_measurement_timescales,
    const DataVector& new_measurement_timescales,
    const int measurements_per_update);

/*!
 * \ingroup ControlSystemGroup
 * \brief Calculate the next expiration time for the MeasurementTimescales.
 *
 * First the function_of_time_expiration_time() is calculated. Then, half of the
 * minimum new measurement timescale is subtracted off from that time and
 * returned. The reason for this is as follows:
 *
 * We update the functions of time one (old) measurement before the expiration
 * time. Based on how dense triggers are set up, which control_system::Trigger
 * is a dense trigger, you calculate the next trigger (measurement) time at the
 * current measurement time. However, at the function of time expiration time we
 * need updated damping timescales from all control systems in order to
 * calculate when the next measurement is going to be (and in turn, the next
 * measurement expiration time). Thus, at the measurement that occurs at the
 * function of time expiration time, our measurement timescales can't be valid
 * and we must wait for updated ones. To achieve this, we set the measurement
 * expiration time *before* the function of time expiration time, but *after*
 * the previous measurement (the update measurement). The factor of one half is
 * just to guarantee we are more than epsilon before the function of time
 * expiration time and more than epsilon after the update measurement.
 */
double measurement_expiration_time(const double time,
                                   const DataVector& old_measurement_timescales,
                                   const DataVector& new_measurement_timescales,
                                   const int measurements_per_update);

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
 * If the control system isn't active then expiration time is
 * `std::numeric_limits<double>::infinity()`.
 *
 * To protect against bad inputs, if the initial expiration time that is
 * calculated is smaller than the initial time step, then the expiration time is
 * simply set to the initial time step. However, the MeasurementTimescales have
 * the same protection so if this does happen, then something is most likely
 * wrong with your initial parameters for the control system.
 */
template <size_t Dim, typename... OptionHolders>
std::unordered_map<std::string, double> initial_expiration_times(
    const double initial_time, const double initial_time_step,
    const int measurements_per_update,
    const std::unique_ptr<::DomainCreator<Dim>>& domain_creator,
    const OptionHolders&... option_holders) {
  std::unordered_map<std::string, double> initial_expiration_times{};

  [[maybe_unused]] const auto gather_initial_expiration_times =
      [&initial_time, &initial_time_step, &measurements_per_update,
       &domain_creator, &initial_expiration_times](const auto& option_holder) {
        const auto& controller = option_holder.controller;
        const std::string& name =
            std::decay_t<decltype(option_holder)>::control_system::name();
        auto tuner = option_holder.tuner;
        Tags::detail::initialize_tuner(make_not_null(&tuner), domain_creator,
                                       initial_time, name);

        const DataVector measurement_timescales =
            calculate_measurement_timescales(controller, tuner,
                                             measurements_per_update);

        const double initial_expiration_time = function_of_time_expiration_time(
            initial_time, DataVector{measurement_timescales.size(), 0.0},
            measurement_timescales, measurements_per_update);
        // Don't have to worry about if functions of time are being overridden
        // because that will be taken care of elsewhere.
        initial_expiration_times[name] =
            option_holder.is_active ? std::max(initial_time + initial_time_step,
                                               initial_expiration_time)
                                    : std::numeric_limits<double>::infinity();
      };

  EXPAND_PACK_LEFT_TO_RIGHT(gather_initial_expiration_times(option_holders));

  return initial_expiration_times;
}
}  // namespace control_system
