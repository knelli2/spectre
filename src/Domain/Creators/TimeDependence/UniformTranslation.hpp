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
#include <utility>
#include <vector>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace FunctionsOfTime {
class FunctionOfTime;
}  // namespace FunctionsOfTime
namespace CoordinateMaps {
namespace TimeDependent {
template <typename Map1, typename Map2, typename Map3>
class ProductOf3Maps;
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
 * \brief A uniform translation in the \f$x-, y-\f$ and \f$z-\f$direction.
 *
 * The coordinates are adjusted according to:
 *
 * \f{align}{
 * x^i \to x^i + f^i(t)
 * \f}
 *
 * where \f$f^i(t)\f$ are the functions of time.
 */
template <size_t MeshDim>
class UniformTranslation final : public TimeDependence<MeshDim> {
 private:
  using TranslationMap =
      domain::CoordinateMaps::TimeDependent::Translation<MeshDim>;

 public:
  using maps_list = tmpl::list<
      domain::CoordinateMap<Frame::Grid, Frame::Inertial, TranslationMap>>;

  static constexpr size_t mesh_dim = MeshDim;

  /// \brief The initial time of the functions of time.
  struct InitialTime {
    using type = double;
    static constexpr Options::String help = {
        "The initial time of the functions of time"};
  };
  /// \brief The \f$x\f$-, \f$y\f$-, and \f$z\f$-velocity.
  struct Velocity {
    using type = std::array<double, MeshDim>;
    static constexpr Options::String help = {"The velocity of the map."};
  };

  using MapForComposition =
      detail::generate_coordinate_map_t<tmpl::list<TranslationMap>>;

  using options = tmpl::list<InitialTime, Velocity>;

  static constexpr Options::String help = {
      "A spatially uniform translation initialized with a constant velocity."};

  UniformTranslation() = default;
  ~UniformTranslation() override = default;
  UniformTranslation(const UniformTranslation&) = delete;
  UniformTranslation(UniformTranslation&&) = default;
  UniformTranslation& operator=(const UniformTranslation&) = delete;
  UniformTranslation& operator=(UniformTranslation&&) = default;

  UniformTranslation(double initial_time,
                     const std::array<double, MeshDim>& velocity);

  auto get_clone() const -> std::unique_ptr<TimeDependence<MeshDim>> override;

  auto block_maps(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, MeshDim>>> override;

  auto functions_of_time(const std::vector<std::pair<std::string, double>>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

  /// Returns the map for each block to be used in a composition of
  /// `TimeDependence`s.
  MapForComposition map_for_composition() const;

 private:
  template <size_t LocalDim>
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend bool operator==(const UniformTranslation<LocalDim>& lhs,
                         const UniformTranslation<LocalDim>& rhs);

  double initial_time_{std::numeric_limits<double>::signaling_NaN()};
  std::array<double, MeshDim> velocity_{};
  mutable std::string function_of_time_name_{};
  mutable double expiration_time_{std::numeric_limits<double>::infinity()};
};

template <size_t Dim>
bool operator==(const UniformTranslation<Dim>& lhs,
                const UniformTranslation<Dim>& rhs);

template <size_t Dim>
bool operator!=(const UniformTranslation<Dim>& lhs,
                const UniformTranslation<Dim>& rhs);
}  // namespace time_dependence
}  // namespace creators
}  // namespace domain
