// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/Creators/TimeDependence/CubicScale.hpp"
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/OptionTags.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformRotationAboutZAxis.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Options/Options.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame
/// \endcond

namespace domain {
namespace creators {
namespace time_dependence {
/*!
 * \brief A TimeDependence that is a composition of a `CubicScale` and a
 * `UniformRotationAboutZAxis`.
 *
 * To create this from options, use something like
 *
 * ```
 * CompositionCubicScaleAndUniformRotationAboutZAxis:
 *   CubicScale:
 *     CubicScale:
 *       InitialTime:
 *       OuterBoundary:
 *       InitialExpansion:
 *       Velocity:
 *       Acceleration:
 *       UseLinearScaling:
 *   UniformRotationAboutZAxis:
 *     UniformRotationAboutZAxis:
 *       InitialTime:
 *       AngularVelocity:
 * ```
 */
template <size_t MeshDim>
class CompositionCubicScaleAndUniformRotationAboutZAxis final
    : public TimeDependence<MeshDim> {
 private:
  using CubicScaleMap =
      domain::CoordinateMaps::TimeDependent::CubicScale<MeshDim>;
  using Rotation = domain::CoordinateMaps::TimeDependent::Rotation<2>;
  using Identity = domain::CoordinateMaps::Identity<1>;
  using RotationMap =
      tmpl::conditional_t<MeshDim == 2, Rotation,
                          domain::CoordinateMaps::TimeDependent::ProductOf2Maps<
                              Rotation, Identity>>;

 public:
  static constexpr size_t mesh_dim = MeshDim;

  using CoordMap = domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                         CubicScaleMap, RotationMap>;

  using maps_list = tmpl::list<CoordMap>;
  static constexpr Options::String help = {
      "A composition of a CubicScale and a UniformRotationAboutZAxis "
      "TimeDependences."};

  using options =
      tmpl::list<OptionTags::TimeDependenceCompositionTag<CubicScale<MeshDim>>,
                 OptionTags::TimeDependenceCompositionTag<
                     UniformRotationAboutZAxis<MeshDim>>>;

  CompositionCubicScaleAndUniformRotationAboutZAxis() = default;
  ~CompositionCubicScaleAndUniformRotationAboutZAxis() override = default;
  CompositionCubicScaleAndUniformRotationAboutZAxis(
      const CompositionCubicScaleAndUniformRotationAboutZAxis&) = default;
  CompositionCubicScaleAndUniformRotationAboutZAxis& operator=(
      const CompositionCubicScaleAndUniformRotationAboutZAxis&) = default;
  CompositionCubicScaleAndUniformRotationAboutZAxis(
      CompositionCubicScaleAndUniformRotationAboutZAxis&&) = default;
  CompositionCubicScaleAndUniformRotationAboutZAxis& operator=(
      CompositionCubicScaleAndUniformRotationAboutZAxis&&) = default;

  explicit CompositionCubicScaleAndUniformRotationAboutZAxis(
      const std::unique_ptr<TimeDependence<MeshDim>>& cubic_scale,
      const std::unique_ptr<TimeDependence<MeshDim>>&
          uniform_rotation_about_z_axis);

  auto get_clone() const -> std::unique_ptr<TimeDependence<MeshDim>> override;

  auto block_maps(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, MeshDim>>> override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

 private:
  std::unique_ptr<TimeDependence<MeshDim>> cubic_scale_;
  std::unique_ptr<TimeDependence<MeshDim>> uniform_rotation_about_z_axis_;

  CoordMap coord_map_;
};
}  // namespace time_dependence
}  // namespace creators
}  // namespace domain
