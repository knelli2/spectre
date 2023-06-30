// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/ControlErrors/Size/AhSpeed.hpp"

#include <cmath>
#include <memory>
#include <sstream>
#include <string>

#include "ControlSystem/ControlErrors/Size/DeltaR.hpp"
#include "Utilities/StdHelpers.hpp"

namespace control_system::size::States {

std::unique_ptr<State> AhSpeed::get_clone() const {
  return std::make_unique<AhSpeed>(*this);
}

std::string AhSpeed::update(const gsl::not_null<Info*> info,
                            const StateUpdateArgs& update_args,
                            const CrossingTimeInfo& crossing_time_info) const {
  const double min_char_speed = update_args.min_char_speed;
  const double min_comoving_char_speed = update_args.min_comoving_char_speed;

  // Note that delta_radius_is_in_danger and char_speed_is_in_danger
  // can be different for different States.
  //
  // The value of 20 causes the control system to panic easily if
  // delta_radius is decreasing quickly.  The value 20 was chosen
  // by trial-and-error in SpEC.
  constexpr double time_tolerance_for_delta_r_in_danger = 20.0;
  const bool delta_radius_is_in_danger =
      crossing_time_info.horizon_will_hit_excision_boundary_first and
      crossing_time_info.t_delta_radius.value_or(
          std::numeric_limits<double>::infinity()) <
          info->damping_time * time_tolerance_for_delta_r_in_danger;

  const bool char_speed_is_in_danger = [&crossing_time_info, &info,
                                        &delta_radius_is_in_danger,
                                        &min_char_speed]() {
    // speed_tolerance_for_char_speed_in_danger is slightly greater
    // than unity so that we don't panic if the actual char speed is
    // large enough compared to the target char speed.  The value 1.1
    // was chosen in SpEC and we haven't had to change it; the behavior
    // should not be sensitive to small changes in this value.
    constexpr double speed_tolerance_for_char_speed_in_danger = 1.1;
    // We don't want to panic unless crossing time is less than
    // the current damping time. The value 0.99 was chosen in SpEC and we
    // haven't had to change it; the behavior should not be sensitive to
    // small changes in this value.
    constexpr double time_tolerance_for_char_speed_in_danger = 0.99;
    return (not delta_radius_is_in_danger) and
           (crossing_time_info.char_speed_will_hit_zero_first and
            crossing_time_info.t_char_speed.value_or(
                std::numeric_limits<double>::infinity()) <
                info->damping_time * time_tolerance_for_char_speed_in_danger and
            min_char_speed < info->target_char_speed *
                                 speed_tolerance_for_char_speed_in_danger);
  }();

  const bool comoving_decreasing_slower_than_char_speeds = not(
      crossing_time_info.t_char_speed.has_value() and
      crossing_time_info.t_comoving_char_speed.has_value() and
      update_args.min_comoving_char_speed > 0.0 and
      update_args.min_comoving_char_speed /
              crossing_time_info.t_comoving_char_speed.value() >
          update_args.min_char_speed / crossing_time_info.t_char_speed.value());
  // The value of 5.0 was chosen by trial and error in SpEC
  constexpr double comoving_char_speed_to_damping_time_limit = 5.0;

  std::stringstream ss{};

  if (char_speed_is_in_danger) {
    ss << "Current state AhSpeed. Char speed in danger. Staying in "
          "AhSpeed.\n";
    if (info->target_char_speed < min_char_speed) {
      // The value of 1.01 below was chosen by trial and error in
      // SpEC.  The value needs to be greater than unity, and something
      // that doesn't drive min_char_speed too quickly.
      constexpr double min_char_speed_increase_factor = 1.01;
      // We are already in state AhSpeed, and we are in danger.
      // But target_char_speed is less than min_char_speed, so we don't want
      // to continue to push min_char_speed downward.  Instead, increase
      // target_char_speed to be above min_char_speed.
      // We don't do this if delta_radius_is_in_danger, because for that
      // case we might need to drive min_char_speed to a smaller value.
      info->target_char_speed = min_char_speed * min_char_speed_increase_factor;
      ss << " Target char speed = " << info->target_char_speed << "\n";
    }
    info->suggested_time_scale = crossing_time_info.t_char_speed;
    ss << " Suggested timescale = " << info->suggested_time_scale;
  } else if (delta_radius_is_in_danger) {
    // The values of target_speed_decrease_factor and
    // time_tolerance_for_delta_r were chosen by trial and error in SpEC.
    constexpr double target_speed_decrease_factor = 0.125;
    constexpr double time_tolerance_for_delta_r =
        0.25 * time_tolerance_for_delta_r_in_danger;
    ss << "Current state AhSpeed. Delta radius in danger.";
    if ((crossing_time_info.t_char_speed.has_value() and
         crossing_time_info.t_delta_radius.value_or(-1.0) >
             info->damping_time * time_tolerance_for_delta_r) or
        update_args.min_comoving_char_speed < 0.0) {
      info->discontinuous_change_has_occurred = true;
      info->target_char_speed = min_char_speed * target_speed_decrease_factor;
      info->suggested_time_scale = std::min(
          info->damping_time, crossing_time_info.t_delta_radius.value_or(
                                  std::numeric_limits<double>::infinity()));

      ss << " Staying in AhSpeed.\n";
      if (update_args.min_comoving_char_speed < 0.0) {
        ss << " Min comoving char speed " << update_args.min_comoving_char_speed
           << " < 0.0\n";
      } else {
        ss << " Char speed X time " << crossing_time_info.t_char_speed.value()
           << " exists and delta radius X time "
           << crossing_time_info.t_delta_radius.value()
           << " > damping time * tolerance "
           << info->damping_time * time_tolerance_for_delta_r << "\n";
      }
      ss << " Target char speed = " << info->target_char_speed << "\n";
      ss << " Suggested timescale = " << info->suggested_time_scale;
    } else {
      info->discontinuous_change_has_occurred = true;
      info->state = std::make_unique<States::DeltaR>();
      info->suggested_time_scale = crossing_time_info.t_delta_radius;
      ss << " Switching to DeltaR.\n";
      ss << " Suggested timescale = " << info->suggested_time_scale;
      // Here is where possible transition to State DeltaRDriftInward will go.
    }
  } else if (update_args.min_comoving_char_speed > 0.0 and
             update_args.min_char_speed > 0.0 and
             (not crossing_time_info.t_comoving_char_speed.has_value() or
              (crossing_time_info.t_comoving_char_speed.value() >
                   comoving_char_speed_to_damping_time_limit *
                       info->damping_time and
               comoving_decreasing_slower_than_char_speeds)) and
             (update_args.min_char_speed >= info->target_char_speed or
              min_comoving_char_speed > min_char_speed)) {
    info->discontinuous_change_has_occurred = true;
    info->state = std::make_unique<States::DeltaR>();

    ss << "Current state AhSpeed. Switching to DeltaR.\n";
    ss << " Min char speed " << update_args.min_char_speed << " > 0\n";
    ss << " Min comoving char speed " << update_args.min_comoving_char_speed
       << " > 0\n";
    if (crossing_time_info.t_comoving_char_speed.has_value()) {
      ss << " Comoving char speed X time "
         << crossing_time_info.t_comoving_char_speed.value()
         << " > damping time * limit "
         << comoving_char_speed_to_damping_time_limit * info->damping_time
         << "\n";
      ss << " Comoving char speeds decreasing slower than char speeds\n";
    }
    if (min_comoving_char_speed > min_char_speed) {
      ss << " Min comoving char speed " << min_comoving_char_speed
         << " > min char speed " << min_char_speed << "\n";
    } else {
      ss << " Min char speed " << update_args.min_char_speed
         << " >= target char speed " << info->target_char_speed << "\n";
    }
    // Here is where possible transition to State DeltaRDriftInward
    // will go.
  } else {
    ss << "Current state AhSpeed. No change necessary. Staying in AhSpeed.";
  }
  // If no 'if's are encountered above, then all the info parameters stay
  // the same as they were.

  return ss.str();
}

double AhSpeed::control_error(
    const Info& info, const ControlErrorArgs& control_error_args) const {
  const double Y00 = sqrt(0.25 / M_PI);
  return (info.target_char_speed - control_error_args.min_char_speed) /
         (Y00 * control_error_args.avg_distorted_normal_dot_unit_coord_vector);
}

PUP::able::PUP_ID AhSpeed::my_PUP_ID = 0;
}  // namespace control_system::size::States
