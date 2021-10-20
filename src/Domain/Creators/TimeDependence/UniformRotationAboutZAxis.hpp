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
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace FunctionsOfTime {
class FunctionOfTime;
}  // namespace FunctionsOfTime
namespace CoordinateMaps {
namespace TimeDependent {
template <typename Map1, typename Map2>
class ProductOf2Maps;
}  // namespace TimeDependent
}  // namespace CoordinateMaps

template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame
/// \endcond

namespace domain {
namespace creators {
namespace time_dependence {
/*!
 * \brief A uniform rotation about the \f$z\f$ axis:
 * \f{eqnarray*}
 * x &\to& x \cos \alpha(t) - y \sin \alpha(t)\text{,} \\
 * y &\to& x \sin \alpha(t) + y \cos \alpha(t)\text{,}
 * \f}
 * where \f$\alpha(t)\f$ is a `domain::FunctionsOfTime::FunctionOfTime`. For 3
 * spatial dimensions, \f$z \to z\f$, and the rotation is implemented as a
 * product of the 2D rotation and an identity map. The rotation is undefined
 * (and therefore unimplemented here) for 1 spatial dimension.
 */
template <size_t MeshDim>
class UniformRotationAboutZAxis final : public TimeDependence<MeshDim> {
  static_assert(
      MeshDim > 1,
      "UniformRotationAboutZAxis<MeshDim> undefined for MeshDim == 1");

 private:
  using Identity = domain::CoordinateMaps::Identity<1>;
  using Rotation = domain::CoordinateMaps::TimeDependent::Rotation<2>;

 public:
  using maps_list = tmpl::list<
      domain::CoordinateMap<Frame::Grid, Frame::Inertial, Rotation>,
      domain::CoordinateMap<
          Frame::Grid, Frame::Inertial,
          CoordinateMaps::TimeDependent::ProductOf2Maps<Rotation, Identity>>>;

  static constexpr size_t mesh_dim = MeshDim;

  /// \brief The initial time of the function of time.
  struct InitialTime {
    using type = double;
    static constexpr Options::String help = {
        "The initial time of the function of time"};
  };
  /// \brief The \f$x\f$-, \f$y\f$-, and \f$z\f$-velocity.
  struct AngularVelocity {
    using type = double;
    static constexpr Options::String help = {
        "The angular velocity of the map."};
  };

  using MapForComposition = detail::generate_coordinate_map_t<
      tmpl::list<tmpl::conditional_t<MeshDim == 2, Rotation,
                                     domain::CoordinateMaps::TimeDependent::
                                         ProductOf2Maps<Rotation, Identity>>>>;

  using options = tmpl::list<InitialTime, AngularVelocity>;

  static constexpr Options::String help = {
      "A spatially uniform rotation about the z axis initialized with a "
      "constant angular velocity."};

  UniformRotationAboutZAxis() = default;
  ~UniformRotationAboutZAxis() override = default;
  UniformRotationAboutZAxis(const UniformRotationAboutZAxis&) = delete;
  UniformRotationAboutZAxis(UniformRotationAboutZAxis&&) = default;
  UniformRotationAboutZAxis& operator=(const UniformRotationAboutZAxis&) =
      delete;
  UniformRotationAboutZAxis& operator=(UniformRotationAboutZAxis&&) = default;

  UniformRotationAboutZAxis(double initial_time, double angular_velocity);

  auto get_clone() const -> std::unique_ptr<TimeDependence<MeshDim>> override;

  auto block_maps(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, MeshDim>>> override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

  /// Returns the map for each block to be used in a composition of
  /// `TimeDependence`s.
  MapForComposition map_for_composition() const;

  static std::string name() {
    return pretty_type::short_name<UniformRotationAboutZAxis>();
  }

 private:
  template <size_t LocalDim>
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend bool operator==(const UniformRotationAboutZAxis<LocalDim>& lhs,
                         const UniformRotationAboutZAxis<LocalDim>& rhs);

  double initial_time_{std::numeric_limits<double>::signaling_NaN()};
  double angular_velocity_{std::numeric_limits<double>::signaling_NaN()};
  std::string function_of_time_name_{"Rotation"};
  mutable double expiration_time_{std::numeric_limits<double>::infinity()};
};

template <size_t Dim>
bool operator==(const UniformRotationAboutZAxis<Dim>& lhs,
                const UniformRotationAboutZAxis<Dim>& rhs);

template <size_t Dim>
bool operator!=(const UniformRotationAboutZAxis<Dim>& lhs,
                const UniformRotationAboutZAxis<Dim>& rhs);
}  // namespace time_dependence
}  // namespace creators
}  // namespace domain
