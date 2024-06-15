// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

namespace domain::creators::time_dependent_options {
/*!
 * \brief Class to be used as an option for initializing expansion map
 * coefficients.
 */
struct ExpansionMapOptions {
  using type = Options::Auto<ExpansionMapOptions, Options::AutoLabel::None>;
  static std::string name() { return "ExpansionMap"; }
  static constexpr Options::String help = {
      "Options for a time-dependent expansion of the coordinates. Specify "
      "'None' to not use this map."};

  struct InitialValues {
    using type =
        std::variant<std::array<double, 3>, FromVolumeFile<names::Expansion>>;
    static constexpr Options::String help = {
        "Initial values for the expansion map, its velocity and "
        "acceleration."};
  };

  struct InitialValuesOuterBoundary {
    using type =
        std::variant<std::array<double, 3>, FromVolumeFile<names::Expansion>>;
    static constexpr Options::String help = {
        "Initial values for the expansion map, its velocity and "
        "acceleration at the outer boundary. Unless you are starting from a "
        "checkpoint or continuing an evolution, this option should likely be "
        "[1.0, 0.0, 0.0] at the start of an evolution"};
  };

  struct DecayTimescaleOuterBoundary {
    using type = double;
    static constexpr Options::String help = {
        "A timescale for how fast the outer boundary expansion approaches its "
        "asymptotic value."};
  };

  struct DecayTimescale {
    using type = Options::Auto<double, Options::AutoLabel::None>;
    static constexpr Options::String help = {
        "If specified, a SettleToConstant function of time will be used for "
        "the expansion map and this number will determine the timescale that "
        "the expansion approaches its asymptotic value. If 'None' is "
        "specified, a PiecewisePolynomial function of time will be used for "
        "the expansion map."};
  };

  struct AsymptoticVelocityOuterBoundary {
    using type = Options::Auto<double, Options::AutoLabel::None>;
    static constexpr Options::String help = {
        "If specified, a FixedSpeedCubic function of time will be used for "
        "the expansion map at the outer boundary and this number will "
        "determine its velocity. If 'None' is "
        "specified, a SettleToConstant function of time will be used for the "
        "expansion map at the outer boundary."};
  };

  using options = tmpl::list<InitialValues, InitialValuesOuterBoundary,
                             DecayTimescaleOuterBoundary, DecayTimescale,
                             AsymptoticVelocityOuterBoundary>;

  ExpansionMapOptions() = default;
  ExpansionMapOptions(
      const std::variant<std::array<double, 3>,
                         FromVolumeFile<names::Expansion>>& expansion_values,
      const std::variant<std::array<double, 3>,
                         FromVolumeFile<names::Expansion>>&
          expansion_outer_boundary_values,
      double decay_timescale_outer_boundary_in,
      std::optional<double> decay_timescale_in,
      std::optional<double> asymptotic_velocity_outer_boundary_in,
      const Options::Context& context = {});

  std::array<DataVector, 3> initial_values{};
  std::array<DataVector, 3> initial_values_outer_boundary{};
  double decay_timescale_outer_boundary{};
  std::optional<double> decay_timescale{};
  std::optional<double> asymptotic_velocity_outer_boundary{};
};
}  // namespace domain::creators::time_dependent_options
