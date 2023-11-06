// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <random>

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Wedge.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Framework/TestHelpers.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/Gsl.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {
namespace {
enum class Surface { Inner, Outer };
void check_spherical_boundary(
    const Wedge& multi_wedge, const std::array<double, 3>& source_coord,
    const std::array<double, 3>& target_coord, const double distorted_radius,
    const double inner_radius, const double outer_radius, const Surface surface,
    const size_t /*zeta_index*/, const double /*sign_zeta_index*/) {
  CAPTURE(source_coord);
  CAPTURE(target_coord);

  double radius = std::numeric_limits<double>::signaling_NaN();
  double expected_transition = std::numeric_limits<double>::signaling_NaN();
  double expected_map_over_radius =
      std::numeric_limits<double>::signaling_NaN();
  double expected_orig_rad_over_rad =
      std::numeric_limits<double>::signaling_NaN();

  if (surface == Surface::Inner) {
    radius = inner_radius;
    expected_transition = 1.0;
    expected_map_over_radius = 1.0 / inner_radius;
    expected_orig_rad_over_rad = inner_radius / magnitude(target_coord);
  } else {
    radius = outer_radius;
    expected_transition = 0.0;
    expected_map_over_radius = 0.0;
    expected_orig_rad_over_rad = 1.0;
  }

  CHECK(multi_wedge(source_coord) == expected_transition);
  CHECK(multi_wedge.map_over_radius(source_coord) ==
        approx(expected_map_over_radius));

  CHECK_ITERABLE_APPROX(
      multi_wedge.gradient(source_coord),
      -source_coord / (radius * (outer_radius - inner_radius)));

  std::optional<double> orig_rad_over_rad =
      multi_wedge.original_radius_over_radius(target_coord, distorted_radius);
  CHECK(orig_rad_over_rad.has_value());
  CHECK(orig_rad_over_rad.value() == approx(expected_orig_rad_over_rad));
}

void check_spherical_boundaries_generic_point(
    const Wedge& multi_wedge, const std::array<double, 3>& source_coord,
    const double distorted_radius, const double inner_radius,
    const double outer_radius, const size_t /*zeta_index*/,
    const double /*sign_zeta_index*/) {
  CAPTURE(source_coord);

  const double radius = magnitude(source_coord);
  CHECK(multi_wedge(source_coord) ==
        approx((outer_radius - radius) / (outer_radius - inner_radius)));
  CHECK(multi_wedge.map_over_radius(source_coord) ==
        approx((outer_radius - radius) /
               ((outer_radius - inner_radius) * radius)));
  CHECK_ITERABLE_APPROX(
      multi_wedge.gradient(source_coord),
      -source_coord / (radius * (outer_radius - inner_radius)));

  const std::array<double, 3> target_coord =
      source_coord * (1.0 - multi_wedge(source_coord) * distorted_radius);
  CAPTURE(target_coord);

  std::optional<double> orig_rad_over_rad =
      multi_wedge.original_radius_over_radius(target_coord, distorted_radius);
  CHECK(orig_rad_over_rad.has_value());
  CHECK(orig_rad_over_rad.value() == approx(radius / magnitude(target_coord)));
}

void check_flat_boundary(const Wedge& multi_wedge,
                         std::array<double, 3> source_coord,
                         std::array<double, 3> target_coord,
                         const double distorted_radius,
                         const double inner_radius, const double outer_radius,
                         const Surface surface, const size_t zeta_index,
                         const double sign_zeta_index) {
  CAPTURE(source_coord);
  CAPTURE(target_coord);

  const double radius = magnitude(source_coord);
  const double inner_distance = inner_radius / sqrt(3.0);
  const double outer_distance = outer_radius / sqrt(3.0);

  CAPTURE(inner_distance);
  CAPTURE(outer_distance);
  double expected_transition = std::numeric_limits<double>::signaling_NaN();
  double expected_map_over_radius =
      std::numeric_limits<double>::signaling_NaN();
  double expected_orig_rad_over_rad =
      std::numeric_limits<double>::signaling_NaN();

  if (surface == Surface::Inner) {
    expected_transition = 1.0;
    expected_map_over_radius = 1.0 / inner_distance;
    expected_orig_rad_over_rad = inner_distance / magnitude(target_coord);
  } else {
    expected_transition = 0.0;
    expected_map_over_radius = 0.0;
    expected_orig_rad_over_rad = 1.0;
  }

  CHECK(multi_wedge(source_coord) == approx(expected_transition));
  CHECK(multi_wedge.map_over_radius(source_coord) ==
        approx(expected_map_over_radius));

  const double one_over_distance_difference =
      1.0 / (outer_distance - inner_distance);
  const auto grad = [&](const double distance) {
    return distance * (square(radius) - square(distance)) / cube(radius);
  };
  // Only for both flat surfaces along z axis at flat surfaces
  const double mag_grad = (grad(outer_distance) - 1.0 -
                           (outer_distance - radius) *
                               (grad(outer_distance) - grad(inner_distance)) *
                               one_over_distance_difference) *
                          one_over_distance_difference;
  std::array<double, 3> expected_gradient{0.0, 0.0, 0.0};
  gsl::at(expected_gradient, zeta_index) = sign_zeta_index * mag_grad;
  CHECK_ITERABLE_APPROX(multi_wedge.gradient(source_coord), expected_gradient);

  std::optional<double> orig_rad_over_rad =
      multi_wedge.original_radius_over_radius(target_coord, distorted_radius);
  CHECK(orig_rad_over_rad.has_value());
  CHECK(orig_rad_over_rad.value() == approx(expected_orig_rad_over_rad));
}

void check_flat_generic_point(const Wedge& /*multi_wedge*/,
                              const std::array<double, 3>& /*source_coord*/,
                              const double /*distorted_radius*/,
                              const double /*inner_radius*/,
                              const double /*outer_radius*/,
                              const size_t /*zeta_index*/,
                              const double /*sign_zeta_index*/) {}

template <typename Generator, typename BoundaryFunction,
          typename GenericFunction>
void test(
    const gsl::not_null<Generator*> /*generator*/,
    const gsl::not_null<
        std::uniform_real_distribution<double>*> /*radial_dist*/,
    const gsl::not_null<std::uniform_real_distribution<double>*> /*polar_dist*/,
    const gsl::not_null<
        std::uniform_real_distribution<double>*> /*azimuthal_dist*/,
    const double inner_radius, const double outer_radius,
    const double inner_sphericity, const double outer_sphericity,
    BoundaryFunction boundary_function, GenericFunction /*generic_function*/,
    const OrientationMap<3>& orientation_map, const size_t zeta_index,
    const double sign_zeta_index) {
  CAPTURE(orientation_map);
  CAPTURE(zeta_index);
  CAPTURE(sign_zeta_index);
  CAPTURE(inner_radius);
  CAPTURE(outer_radius);
  CAPTURE(inner_sphericity);
  CAPTURE(outer_sphericity);

  Wedge multi_wedge{inner_radius, outer_radius, inner_sphericity,
                    outer_sphericity, orientation_map};

  // The distortion can't be 1 or larger, otherwise the coordinates will be
  // mapped through the center of the grid. The original_radius_over_radius
  // function starts going bad as distorted_radius approaches 1.
  const double distorted_radius = 0.8;

  std::array<double, 3> source_coord_aligned{0.0, 0.0, inner_radius};
  // The orientation map takes a coordinate that's aligned in the +zeta
  // direction and maps it to it's expected orientation. Thus we want to do this
  // for each orientation so the source coord is properly aligned
  std::array<double, 3> source_coord = source_coord_aligned;
  const auto set_source_coord =
      [&](const std::array<double, 3>& new_source_coord_aligned) {
        source_coord_aligned = new_source_coord_aligned;
        source_coord = discrete_rotation(orientation_map, source_coord_aligned);
      };

  {
    INFO("z-axis inner radius");
    set_source_coord(inner_sphericity == 1.0 ? source_coord
                                             : source_coord / sqrt(3.0));
    CAPTURE(source_coord_aligned);
    // Assume center is 0 and sum of Ylm is distorted radius above. Transition
    // function at inner radius is 1.
    const std::array<double, 3> target_coord =
        source_coord * (1.0 - distorted_radius);

    boundary_function(multi_wedge, source_coord, target_coord, distorted_radius,
                      inner_radius, outer_radius, Surface::Inner, zeta_index,
                      sign_zeta_index);
  }

  // {
  //   INFO("z-axis outer radius");
  //   set_source_coord(
  //       {0.0, 0.0,
  //        outer_sphericity == 1.0 ? outer_radius : outer_radius / sqrt(3.0)});
  //   CAPTURE(source_coord_aligned);

  //   // Target coord is source coord
  //   boundary_function(multi_wedge, source_coord, source_coord,
  //   distorted_radius,
  //                     inner_radius, outer_radius, Surface::Outer, zeta_index,
  //                     sign_zeta_index);
  // }

  // {
  //   INFO("z-axis radius between inner and outer");
  //   const double radius = (*radial_dist)(*generator);
  //   set_source_coord({0.0, 0.0, radius});
  //   CAPTURE(source_coord_aligned);

  //   generic_function(multi_wedge, source_coord, distorted_radius,
  //   inner_radius,
  //                    outer_radius, zeta_index, sign_zeta_index);
  // }

  // {
  //   INFO("Corners inner surface");
  //   const double length = inner_radius / sqrt(3.0);
  //   for (const auto& [x, y] : cartesian_product(std::array{length, -length},
  //                                               std::array{length, -length}))
  //                                               {
  //     set_source_coord({x, y, length});
  //     CAPTURE(source_coord_aligned);
  //     // Assume center is 0 and sum of Ylm is distorted radius above.
  //     Transition
  //     // function at inner radius is 1.
  //     const std::array<double, 3> target_coord =
  //         source_coord * (1.0 - distorted_radius);

  //     boundary_function(multi_wedge, source_coord, target_coord,
  //                       distorted_radius, inner_radius, outer_radius,
  //                       Surface::Inner, zeta_index, sign_zeta_index);
  //   }
  // }

  // {
  //   INFO("Corners outer surface");
  //   const double length = outer_radius / sqrt(3.0);
  //   for (const auto& [x, y] : cartesian_product(std::array{length, -length},
  //                                               std::array{length, -length}))
  //                                               {
  //     set_source_coord({x, y, length});
  //     CAPTURE(source_coord_aligned);

  //     // Target coord is source coord
  //     boundary_function(multi_wedge, source_coord, source_coord,
  //                       distorted_radius, inner_radius, outer_radius,
  //                       Surface::Outer, zeta_index, sign_zeta_index);
  //   }
  // }

  // {
  //   INFO("Corners between inner and outer");
  //   const double radius = (*radial_dist)(*generator);
  //   const double length = radius / sqrt(3.0);
  //   for (const auto& [x, y] : cartesian_product(std::array{length, -length},
  //                                               std::array{length, -length}))
  //                                               {
  //     set_source_coord({x, y, length});
  //     CAPTURE(source_coord_aligned);

  //     generic_function(multi_wedge, source_coord, distorted_radius,
  //                      inner_radius, outer_radius, zeta_index,
  //                      sign_zeta_index);
  //   }
  // }

  // {
  //   INFO("Random angle, random radius");
  //   for (size_t i = 0; i < 10; i++) {
  //     const double radius = (*radial_dist)(*generator);
  //     const double theta = (*polar_dist)(*generator);
  //     const double phi = (*azimuthal_dist)(*generator);
  //     set_source_coord({radius * sin(theta) * cos(phi),
  //                       radius * sin(theta) * sin(phi), radius *
  //                       cos(theta)});
  //     CAPTURE(source_coord_aligned);

  //     generic_function(multi_wedge, source_coord, distorted_radius,
  //                      inner_radius, outer_radius, zeta_index,
  //                      sign_zeta_index);
  //   }
  // }
}

SPECTRE_TEST_CASE("Unit.Domain.CoordinateMaps.Shape.Wedge", "[Domain][Unit]") {
  MAKE_GENERATOR(generator);
  const double inner_radius = 1.5;
  const double outer_radius = 8.0;
  std::uniform_real_distribution<double> radial_dist{inner_radius,
                                                     outer_radius};
  // Wedge only opens pi/2 in both directions
  std::uniform_real_distribution<double> polar_dist{0.0, M_PI_4};
  std::uniform_real_distribution<double> azimuthal_dist{0.0, 2.0 * M_PI};

  const auto orientations = orientations_for_sphere_wrappings();
  // This is based off ordering of orientations. From DomainHelpers.cpp it goes
  // +z, -z, +y, -y, +x, -x
  size_t zeta_index = 2;
  double sign_zeta_index = 1.0;
  for (size_t i = 0; i < orientations.size(); i++) {
    test(make_not_null(&generator), make_not_null(&radial_dist),
         make_not_null(&polar_dist), make_not_null(&azimuthal_dist),
         inner_radius, outer_radius, 1.0, 1.0, check_spherical_boundary,
         check_spherical_boundaries_generic_point, gsl::at(orientations, i),
         zeta_index, sign_zeta_index);
    test(make_not_null(&generator), make_not_null(&radial_dist),
         make_not_null(&polar_dist), make_not_null(&azimuthal_dist),
         inner_radius, outer_radius, 0.0, 0.0, check_flat_boundary,
         check_flat_generic_point, gsl::at(orientations, i), zeta_index,
         sign_zeta_index);
    zeta_index -= i % 2;
    sign_zeta_index *= -1.0;
  }
  (void)check_flat_boundary;
}
}  // namespace
}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
