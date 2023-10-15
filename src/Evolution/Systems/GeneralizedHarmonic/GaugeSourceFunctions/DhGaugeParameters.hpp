// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <string>

#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace gh::gauges {
/*!
 * \brief A struct holding the parameters for initializing damped harmonic gauge
 *
 * If UseRollon is true, the gauge transitions ("rolls on") to damped
 * harmonic from the initial data gauge; otherwise, the gauge begins
 * immediately in damped harmonic gauge.
 */
template <bool UseRollon>
struct DhGaugeParameters;

/// \cond
template <>
struct DhGaugeParameters<true> {
  double rollon_start;
  double rollon_window;
  double spatial_decay_width;
  std::array<double, 3> amplitudes;
  std::array<int, 3> exponents;

  inline const static std::string help{
      "A struct holding the parameters for initializing damped harmonic "
      "gauge, including a roll-on from the initial gauge."};

  /// The rollon start time
  struct RollOnStartTime {
    using type = double;
    inline const static std::string help{
        "Simulation time to start rolling on the damped harmonic gauge"};
  };

  /// The width of the Gaussian for the gauge rollon
  struct RollOnTimeWindow {
    using type = double;
    inline const static std::string help{
        "The width of the Gaussian that controls how quickly the gauge is "
        "rolled on."};
  };

  /// The width of the Gaussian for the spatial decay of the damped harmonic
  /// gauge.
  struct SpatialDecayWidth {
    using type = double;
    inline const static std::string help{
        "Spatial width of weight function used in the damped harmonic "
        "gauge."};
  };

  /// The amplitudes for the L1, L2, and S terms, respectively, for the damped
  /// harmonic gauge.
  struct Amplitudes {
    using type = std::array<double, 3>;
    inline const static std::string help{
        "Amplitudes [AL1, AL2, AS] for the damped harmonic gauge."};
  };

  /// The exponents for the L1, L2, and S terms, respectively, for the damped
  /// harmonic gauge.
  struct Exponents {
    using type = std::array<int, 3>;
    inline const static std::string help{
        "Exponents [eL1, eL2, eS] for the damped harmonic gauge."};
  };

  using options = tmpl::list<RollOnStartTime, RollOnTimeWindow,
                             SpatialDecayWidth, Amplitudes, Exponents>;

  DhGaugeParameters(double start, double window, double width,
                    const std::array<double, 3>& amps,
                    const std::array<int, 3>& exps);

  DhGaugeParameters() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};

template <>
struct DhGaugeParameters<false> {
  double spatial_decay_width;
  std::array<double, 3> amplitudes;
  std::array<int, 3> exponents;

  inline const static std::string help{
      "A struct holding the parameters for initializing damped harmonic "
      "gauge with no roll-on from the initial gauge."};

  /// The width of the Gaussian for the spatial decay of the damped harmonic
  /// gauge.
  struct SpatialDecayWidth {
    using type = double;
    inline const static std::string help{
        "Spatial width of weight function used in the damped harmonic "
        "gauge."};
  };

  /// The amplitudes for the L1, L2, and S terms, respectively, for the damped
  /// harmonic gauge.
  struct Amplitudes {
    using type = std::array<double, 3>;
    inline const static std::string help{
        "Amplitudes [AL1, AL2, AS] for the damped harmonic gauge."};
  };

  /// The exponents for the L1, L2, and S terms, respectively, for the damped
  /// harmonic gauge.
  struct Exponents {
    using type = std::array<int, 3>;
    inline const static std::string help{
        "Exponents [eL1, eL2, eS] for the damped harmonic gauge."};
  };

  using options = tmpl::list<SpatialDecayWidth, Amplitudes, Exponents>;

  DhGaugeParameters(double width, const std::array<double, 3>& amps,
                    const std::array<int, 3>& exps);

  DhGaugeParameters() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);
};
/// \endcond
}  // namespace gh::gauges
