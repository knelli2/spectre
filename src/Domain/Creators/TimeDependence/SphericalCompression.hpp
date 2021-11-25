// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/TimeDependent/SphericalCompression.hpp"
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace FunctionsOfTime {
class FunctionOfTime;
}  // namespace FunctionsOfTime
namespace CoordinateMaps::TimeDependent {
template <bool InertiorMap>
class SphericalCompression;
}  // namespace CoordinateMaps::TimeDependent
template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame
/// \endcond

namespace domain::creators::time_dependence {
/*!
 * \brief A spherical compression about some center, as given by
 * domain::CoordinateMaps::TimeDependent::SphericalCompression<false>.
 *
 * \details This TimeDependence is suitable for use on a spherical shell,
 * where MinRadius and MaxRadius are the inner and outer radii of the shell,
 * respectively.
 */
class SphericalCompression final : public TimeDependence<3> {
 private:
  using SphericalCompressionMap =
      domain::CoordinateMaps::TimeDependent::SphericalCompression<false>;

 public:
  using maps_list =
      tmpl::list<domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                       SphericalCompressionMap>>;

  static constexpr size_t mesh_dim = 3;

  /// \brief The initial time of the function of time.
  struct InitialTime {
    using type = double;
    static constexpr Options::String help = {
        "The initial time of the function of time"};
  };
  /// \brief Minimum radius for the SphericalCompression map
  struct MinRadius {
    using type = double;
    static constexpr Options::String help = {
        "Min radius for SphericalCompression map."};
  };
  /// \brief Maximum radius for the SphericalCompression map
  struct MaxRadius {
    using type = double;
    static constexpr Options::String help = {
        "Max radius for SphericalCompression map."};
  };
  /// \brief Center for the SphericalCompression map
  struct Center {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "Center for the SphericalCompression map."};
  };
  /// \brief Initial value for function of time for the spherical compression
  struct InitialValue {
    using type = double;
    static constexpr Options::String help = {
        "Spherical compression value at initial time."};
  };
  /// \brief Initial radial velocity for the function of time for the spherical
  /// compression
  struct InitialVelocity {
    using type = double;
    static constexpr Options::String help = {
        "Spherical compression initial radial velocity."};
  };
  /// \brief Initial radial acceleration for the function of time for the
  /// spherical compression
  struct InitialAcceleration {
    using type = double;
    static constexpr Options::String help = {
        "Spherical compression initial radial acceleration."};
  };

  using MapForComposition =
      detail::generate_coordinate_map_t<tmpl::list<SphericalCompressionMap>>;

  using options =
      tmpl::list<InitialTime, MinRadius, MaxRadius, Center, InitialValue,
                 InitialVelocity, InitialAcceleration>;

  static constexpr Options::String help = {"A spherical compression."};

  SphericalCompression() = default;
  ~SphericalCompression() override = default;
  SphericalCompression(const SphericalCompression&) = delete;
  SphericalCompression(SphericalCompression&&) = default;
  SphericalCompression& operator=(const SphericalCompression&) = delete;
  SphericalCompression& operator=(SphericalCompression&&) = default;

  SphericalCompression(double initial_time, double min_radius,
                       double max_radius, std::array<double, 3> center,
                       double initial_value, double initial_velocity,
                       double initial_acceleration,
                       const Options::Context& context = {});

  auto get_clone() const -> std::unique_ptr<TimeDependence<mesh_dim>> override;

  auto block_maps(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, mesh_dim>>> override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

  /// Returns the map for each block to be used in a composition of
  /// `TimeDependence`s.
  MapForComposition map_for_composition() const;

 private:
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend bool operator==(const SphericalCompression& lhs,
                         const SphericalCompression& rhs);

  double initial_time_{std::numeric_limits<double>::signaling_NaN()};
  double min_radius_{std::numeric_limits<double>::signaling_NaN()};
  double max_radius_{std::numeric_limits<double>::signaling_NaN()};
  std::array<double, 3> center_{};
  double initial_value_{std::numeric_limits<double>::signaling_NaN()};
  double initial_velocity_{std::numeric_limits<double>::signaling_NaN()};
  double initial_acceleration_{std::numeric_limits<double>::signaling_NaN()};
  std::string function_of_time_name_{"Size"};
};

bool operator!=(const SphericalCompression& lhs,
                const SphericalCompression& rhs);
}  // namespace domain::creators::time_dependence
