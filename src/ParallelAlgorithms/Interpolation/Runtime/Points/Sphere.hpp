// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Structure/Element.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
struct FunctionOfTime;
template <size_t Dim>
struct Domain;
template <size_t Dim>
struct Element;
}  // namespace domain::Tags
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp2::points {
/*!
 * \brief A series of concentric spherical surfaces.
 *
 * \details The parameter `LMax` sets the number of collocation points on
 * each spherical surface equal to `(l_max + 1) * (2 * l_max + 1)`. The
 * parameter `AngularOrdering` encodes the collocation ordering. For example,
 * the apparent horizon finder relies on spherepack routines that require
 * `Strahlkorper` for `AngularOrdering`, and using these surfaces for a CCE
 * worldtube requires `Cce` for `AngularOrdering`.
 */
struct Sphere : public tt::ConformsTo<protocols::Points> {
 private:
  using BlockCoords = std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, 3, ::Frame::BlockLogical>>>>;

 public:
  using tags_on_target = tmpl::list<>;
  using points_volume_compute_tags = tmpl::list<>;

  struct LMax {
    using type = size_t;
    static constexpr Options::String help = {
        "The number of collocation points on each sphere will be equal to "
        "`(l_max + 1) * (2 * l_max + 1)`"};
  };
  struct Center {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {"Center of every sphere"};
  };
  struct Radius {
    using type = std::variant<double, std::vector<double>>;
    static constexpr Options::String help = {"Radius of the sphere(s)"};
  };
  struct AngularOrdering {
    using type = intrp2::AngularOrdering;
    static constexpr Options::String help = {
        "Chooses theta,phi ordering in 2d array"};
  };
  using options = tmpl::list<LMax, Center, Radius, AngularOrdering>;
  static constexpr Options::String help = {
      "An arbitrary number of spherical surfaces."};

  Sphere(size_t l_max, const std::array<double, 3>& center,
         const typename Radius::type& input_radii,
         intrp2::AngularOrdering angular_ordering,
         const Options::Context& context = {});

  Sphere() = default;

  using argument_tags =
      tmpl::list<domain::Tags::Element<3>,
                 domain::Tags::Coordinates<3, Frame::Grid>,
                 domain::Tags::Coordinates<3, Frame::Distorted>,
                 domain::Tags::Coordinates<3, Frame::Inertial>>;

  template <typename Metavariables>
  std::optional<BlockCoords> operator()(
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& frame, const Element<3>& element,
      const tnsr::I<DataVector, 3, Frame::Grid> grid_coordinates,
      const tnsr::I<DataVector, 3, Frame::Distorted> distorted_coordinates,
      const tnsr::I<DataVector, 3, Frame::Inertial> inertial_coordinates)
      const {
    BlockCoords block_logical_coords{};
    // Too big, but we'll shrink later
    block_logical_coords.reserve(get<0>(points_).size());

    // Extremum for x, y, z, r
    std::array<std::pair<double, double>, 4> min_max_coordinates{};

    const tnsr::I<DataVector, 3, Frame::NoFrame> coordinates{};
    for (size_t i = 0; i < 3; i++) {
      if (frame == "Grid") {
        make_const_view(make_not_null(&coordinates.get(i)),
                        grid_coordinates.get(i), 0,
                        grid_coordinates.get(i).size());
      } else if (frame == "Distorted") {
        make_const_view(make_not_null(&coordinates.get(i)),
                        distorted_coordinates.get(i), 0,
                        distorted_coordinates.get(i).size());
      } else {
        make_const_view(make_not_null(&coordinates.get(i)),
                        inertial_coordinates.get(i), 0,
                        inertial_coordinates.get(i).size());
      }
    }

    // Calculate r^2 from center of sphere because sqrt is expensive
    DataVector radii_squared{get<0>(coordinates).size(), 0.0};
    for (size_t i = 0; i < 3; i++) {
      radii_squared += square(coordinates.get(i) - gsl::at(center_, i));
    }

    // Compute min and max
    {
      const auto [min, max] = alg::minmax_element(radii_squared);
      min_max_coordinates[3].first = *min;
      min_max_coordinates[3].second = *max;
    }

    const size_t number_of_angular_points = (l_max_ + 1) * (2 * l_max_ + 1);
    // first size_t = position of first radius in bounds
    // second size_t = total bounds to use/check
    std::optional<std::pair<size_t, size_t>> offset_and_num_points{};

    // Have a very small buffer just in case of roundoff
    double epsilon =
        (min_max_coordinates[3].second - min_max_coordinates[3].first) *
        std::numeric_limits<double>::epsilon() * 100.0;
    size_t offset_index = 0;
    // Check if any radii of the target are within the radii of our element
    for (double radius : radii_) {
      const double square_radius = square(radius);
      if (square_radius >= (gsl::at(min_max_coordinates, 3).first - epsilon) and
          square_radius <= (gsl::at(min_max_coordinates, 3).second + epsilon)) {
        if (offset_and_num_points.has_value()) {
          offset_and_num_points->second += number_of_angular_points;
        } else {
          offset_and_num_points =
              std::make_pair(offset_index * number_of_angular_points,
                             number_of_angular_points);
        }
      }
      offset_index++;
    }

    // If no radii pass through this element, there's nothing to do so return
    if (not offset_and_num_points.has_value()) {
      return std::nullopt;
    }

    // Get the x,y,z bounds
    for (size_t i = 0; i < 3; i++) {
      const auto [min, max] = alg::minmax_element(coordinates.get(i));
      gsl::at(min_max_coordinates, i).first = *min;
      gsl::at(min_max_coordinates, i).second = *max;
    }

    const tnsr::I<DataVector, 3, Frame::NoFrame> target_points_to_check{};
    // Use the offset and number of points to create a view. We assume that if
    // there are multiple radii in this element, that they are successive
    // radii. If this wasn't true, that'd be a really weird topology.
    for (size_t i = 0; i < 3; i++) {
      make_const_view(make_not_null(&target_points_to_check.get(i)),
                      points_.get(i), offset_and_num_points->first,
                      offset_and_num_points->second);
    }

    // To break out of inner loop and skip the point
    bool skip_point = false;
    tnsr::I<double, 3, Frame::NoFrame> sphere_coords_to_map{};

    const Block<3>& block = Parallel::get<domain::Tags::Domain<3>>(cache)
                                .blocks()[element.id().block_id()];

    // Now for every radius in this element, we check if their points are
    // within the x,y,z bounds of the element. If they are, map the point to
    // the block logical frame and add it to the vector f all block logical
    // coordinates.
    for (size_t index = 0; index < get<0>(target_points_to_check).size();
         index++) {
      skip_point = false;
      for (size_t i = 0; i < 3; i++) {
        const double coord = target_points_to_check.get(i)[index];
        epsilon = (gsl::at(min_max_coordinates, i).second -
                   gsl::at(min_max_coordinates, i).first) *
                  std::numeric_limits<double>::epsilon() * 100.0;
        // If a point is outside any of the bounding box, skip it
        if (coord < (gsl::at(min_max_coordinates, i).first - epsilon) or
            coord > (gsl::at(min_max_coordinates, i).second + epsilon)) {
          skip_point = true;
          break;
        }

        sphere_coords_to_map.get(i) = coord;
      }

      if (skip_point) {
        continue;
      }

      std::optional<tnsr::I<double, 3, ::Frame::BlockLogical>>
          block_coords_of_target_point{};

      if constexpr (Parallel::is_in_global_cache<
                        Metavariables, domain::Tags::FunctionsOfTime>) {
        const auto& functions_of_time =
            Parallel::get<domain::Tags::FunctionsOfTime>(cache);

        block_coords_of_target_point =
            block_logical_coordinates_single_point_in_frame(
                sphere_coords_to_map, frame, block, time, functions_of_time);
      } else {
        block_coords_of_target_point =
            block_logical_coordinates_single_point_in_frame(
                sphere_coords_to_map, frame, block);
      }

      if (block_coords_of_target_point.has_value()) {
        // Get index into vector of all grid points of the target. This is
        // just the offset + index
        block_logical_coords.emplace_back(
            domain::BlockId(element.id().block_id()),
            std::move(block_coords_of_target_point.value()));
      }
    }

    if (block_logical_coords.empty()) {
      return std::nullopt;
    } else {
      return {block_logical_coords};
    }
  }

  /// @{
  /*!
   * \brief Methods specific to the Sphere target that return the input
   * parameters.
   */
  size_t l_max() const;
  const std::array<double, 3>& center() const;
  const std::set<double>& radii() const;
  intrp2::AngularOrdering angular_ordering() const;
  /// @}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  const tnsr::I<DataVector, 3, Frame::NoFrame>& target_points_no_frame() const;

 private:
  size_t l_max_{};
  std::array<double, 3> center_{};
  std::set<double> radii_{};
  intrp2::AngularOrdering angular_ordering_{};
  tnsr::I<DataVector, 3, Frame::NoFrame> points_{};

  friend bool operator==(const Sphere& lhs, const Sphere& rhs);
};

template <size_t Dim>
bool operator!=(const Sphere& lhs, const Sphere& rhs);
}  // namespace intrp2::points
