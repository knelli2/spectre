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
#include "Domain/CoordinateMaps/TimeDependent/RadialTranslation.hpp"
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::FunctionsOfTime {
class FunctionOfTime;
}  // namespace domain::FunctionsOfTime
namespace domain::CoordinateMaps::TimeDependent {
template <typename Map1, typename Map2, typename Map3>
class ProductOf3Maps;
template <typename Map1, typename Map2>
class ProductOf2Maps;
}  // namespace domain::CoordinateMaps::TimeDependent
namespace domain {
template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain

namespace Frame {
struct Distorted;
struct Grid;
struct Inertial;
}  // namespace Frame
/// \endcond

namespace domain::creators::time_dependence {
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
 *
 * \p Index is used to distinguish multiple `RadialTranslation`s from each
 * other in CompositionRadialTranslation.
 *
 * See the documentation for the constructors below: one constructor
 * takes two velocities, which correspond to two translations: one from
 * Frame::Grid to Frame::Distorted, and the other from Frame::Distorted to
 * Frame::Inertial.
 */
template <size_t MeshDim>
class RadialTranslation final : public TimeDependence<MeshDim> {
 private:
  using RadialTranslationMap =
      domain::CoordinateMaps::TimeDependent::RadialTranslation<MeshDim>;

 public:
  using maps_list =
      tmpl::list<domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                       RadialTranslationMap>>;

  static constexpr size_t mesh_dim = MeshDim;

  /// \brief The initial time of the functions of time.
  struct InitialTime {
    using type = double;
    static constexpr Options::String help = {
        "The initial time of the functions of time"};
  };
  struct InnerMapInitialValues {
    using type = std::vector<double>;
    static constexpr Options::String help = {
        "Initial values of the inner function of time. Can also specify "
        "velocity and acceleration."};
    static size_t lower_bound_on_size() { return 1; }
  };
  struct OuterMapOptions {
    using type = Options::Auto<OuterMapOptions, Options::AutoLabel::None>;
    static constexpr Options::String help = {
        "Options for have two functions of time controlling the radial "
        "translation. Specify 'None' to only have a single rigid radial "
        "translation."};

    struct InitialValues {
      using type = std::vector<double>;
      static constexpr Options::String help = {
          "Initial values of the of the outer function of time. Can also "
          "specify velocity and acceleration."};
      static size_t lower_bound_on_size() { return 1; }
    };
    struct InnerOuterRadius {
      using type = std::array<double, 2>;
      static constexpr Options::String help = {
          "Inner and outer radius of the region of radial translation in grid "
          "coordinates."};
    };

    using options = tmpl::list<InitialValues, InnerOuterRadius>;

    OuterMapOptions() = default;
    OuterMapOptions(std::vector<double> initial_values_in,
                    const std::array<double, 2>& radii)
        : initial_values(std::move(initial_values_in)),
          inner_outer_radius(radii) {}

    std::vector<double> initial_values{};
    std::array<double, 2> inner_outer_radius{};
  };

  using GridToInertialMap =
      detail::generate_coordinate_map_t<Frame::Grid, Frame::Inertial,
                                        tmpl::list<RadialTranslationMap>>;

  using options =
      tmpl::list<InitialTime, InnerMapInitialValues, OuterMapOptions>;

  static constexpr Options::String help = {"A radial translation."};

  RadialTranslation() = default;
  ~RadialTranslation() override = default;
  RadialTranslation(const RadialTranslation&) = delete;
  RadialTranslation(RadialTranslation&&) = default;
  RadialTranslation& operator=(const RadialTranslation&) = delete;
  RadialTranslation& operator=(RadialTranslation&&) = default;

  RadialTranslation(double initial_time,
                    std::vector<double> inner_map_initial_values,
                    const std::optional<OuterMapOptions>& outer_map_options,
                    const Options::Context& context = {});

  auto get_clone() const -> std::unique_ptr<TimeDependence<MeshDim>> override;

  auto block_maps_grid_to_inertial(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, MeshDim>>> override;

  auto block_maps_grid_to_distorted(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Distorted, MeshDim>>> override;

  auto block_maps_distorted_to_inertial(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Distorted, Frame::Inertial, MeshDim>>> override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

 private:
  template <size_t LocalDim>
  // NOLINTNEXTLINE(readability-redundant-declaration)
  friend bool operator==(const RadialTranslation<LocalDim>& lhs,
                         const RadialTranslation<LocalDim>& rhs);

  GridToInertialMap grid_to_inertial_map() const;

  double initial_time_{std::numeric_limits<double>::signaling_NaN()};
  std::vector<double> inner_initial_values_{};
  std::optional<OuterMapOptions> outer_map_options_{};
  inline static const std::string inner_function_of_time_name_{
      "InnerRadialTranslation"};
  inline static const std::string outer_function_of_time_name_{
      "OuterRadialTranslation"};
};

template <size_t Dim>
bool operator==(const RadialTranslation<Dim>& lhs,
                const RadialTranslation<Dim>& rhs);

template <size_t Dim>
bool operator!=(const RadialTranslation<Dim>& lhs,
                const RadialTranslation<Dim>& rhs);
}  // namespace domain::creators::time_dependence
