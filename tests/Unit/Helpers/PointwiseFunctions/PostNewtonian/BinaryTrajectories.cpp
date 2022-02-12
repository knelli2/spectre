// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/PointwiseFunctions/PostNewtonian/BinaryTrajectories.hpp"

#include <array>
#include <cmath>
#include <utility>

#include "Utilities/ConstantExpressions.hpp"

BinaryTrajectories::BinaryTrajectories(double initial_separation,
                                       const std::array<double, 3>& velocity,
                                       bool newtonian)
    : initial_separation_fourth_power_{square(square(initial_separation))},
      velocity_(velocity),
      newtonian_(newtonian) {}

double BinaryTrajectories::separation(const double time) const {
  const double pn_correction_term = newtonian_ ? 0.0 : 12.8 * time;
  return pow(initial_separation_fourth_power_ - pn_correction_term, 0.25);
}

double BinaryTrajectories::orbital_frequency(const double time) const {
  return pow(separation(time), -1.5);
}

double BinaryTrajectories::omega(const double time) const {
  // This is d/dt(orbital_frequency) if we are using PN, but 0 if it's newtonian
  const double pn_correction_term =
      newtonian_
          ? 0.0
          : 4.8 * pow(initial_separation_fourth_power_ - 12.8 * time, -1.375);
  // const double pn_correction_term = 0.0;
  return orbital_frequency(time) + pn_correction_term * time;
}

std::pair<std::array<double, 3>, std::array<double, 3>>
BinaryTrajectories::positions(const double time) const {
  const double sep = separation(time);
  const double omega = orbital_frequency(time);
  double xB = 0.5 * sep * cos(omega * time);
  double yB = 0.5 * sep * sin(omega * time);
  double xA = -xB;
  double yA = -yB;
  xB += velocity_[0] * time;
  xA += velocity_[0] * time;
  yB += velocity_[1] * time;
  yA += velocity_[1] * time;
  const double zB = velocity_[2] * time;
  const double zA = velocity_[2] * time;
  return {{xA, yA, zA}, {xB, yB, zB}};
}

std::pair<std::array<double, 3>, std::array<double, 3>>
BinaryTrajectories::positions_no_expansion(const double time) const {
  // Separation stays constant while omega follows PN (or newtonian) values
  const double sep = pow(initial_separation_fourth_power_, 0.25);
  const double omega = orbital_frequency(time);
  double xB = 0.5 * sep * cos(omega * time);
  double yB = 0.5 * sep * sin(omega * time);
  double xA = -xB;
  double yA = -yB;
  xB += velocity_[0] * time;
  xA += velocity_[0] * time;
  yB += velocity_[1] * time;
  yA += velocity_[1] * time;
  const double zB = velocity_[2] * time;
  const double zA = velocity_[2] * time;
  return {{xA, yA, zA}, {xB, yB, zB}};
}
