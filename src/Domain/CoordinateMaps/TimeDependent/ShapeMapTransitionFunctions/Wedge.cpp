// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Wedge.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <pup.h>

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "NumericalAlgorithms/RootFinding/QuadraticEquation.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"

#include "Domain/Structure/Direction.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/ErrorHandling/CaptureForError.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {
template <typename T>
T Wedge::Surface::distance(const std::array<T, 3>& coords) const {
  // Short circuit if it's a sphere. Then the distance is trivially the radius
  // of this surface
  if (sphericity == 1.0) {
    return make_with_value<T>(coords[0], radius);
  }

  // D = R * [ (1 - s) / (sqrt(3) * cos(theta)) + s]
  // cos(theta) = z / r
  return radius *
         ((1.0 - sphericity) / (coords[2] * sqrt(3.0)) * magnitude(coords) +
          sphericity);
}

void Wedge::Surface::pup(PUP::er& p) {
  p | radius;
  p | sphericity;
}

bool Wedge::Surface::operator==(const Wedge::Surface& other) const {
  return radius == other.radius and sphericity == other.sphericity;
}

bool Wedge::Surface::operator!=(const Wedge::Surface& other) const {
  return not(*this == other);
}

Wedge::Wedge(double overall_inner_radius, double overall_outer_radius,
             double overall_inner_sphericity, double overall_outer_sphericity,
             double this_wedge_inner_radius, double this_wedge_outer_radius,
             double this_wedge_inner_sphericity,
             double this_wedge_outer_sphericity,
             OrientationMap<3> orientation_map)
    : overall_inner_surface_(
          Surface{overall_inner_radius, overall_inner_sphericity}),
      overall_outer_surface_(
          Surface{overall_outer_radius, overall_outer_sphericity}),
      this_wedge_inner_surface_(
          Surface{this_wedge_inner_radius, this_wedge_inner_sphericity}),
      this_wedge_outer_surface_(
          Surface{this_wedge_outer_radius, this_wedge_outer_sphericity}),
      orientation_map_(std::move(orientation_map)) {
  direction_ = orientation_map_.mapped_directions()[1];
  if (direction_.axis() == Direction<3>::Axis::Eta) {
    direction_ = Direction<3>{Direction<3>::Axis::Zeta, direction_.side()};
  } else if (direction_.axis() == Direction<3>::Axis::Zeta) {
    direction_ = Direction<3>{Direction<3>::Axis::Eta, direction_.side()};
  }
  Parallel::printf(
      "Constructing Wedge transition with\n"
      " overall r_in:  %.16f\n"
      " overall r_out: %.16f\n"
      " overall s_in:  %f\n"
      " overall s_out: %f\n"
      " this wedge's r_in:  %.16f\n"
      " this wedge's r_out: %.16f\n"
      " this wedge's s_in:  %f\n"
      " this wedge's s_out: %f\n"
      " orientation map: %s\n"
      " direction: %s\n",
      overall_inner_radius, overall_outer_radius, overall_inner_sphericity,
      overall_outer_sphericity, this_wedge_inner_radius,
      this_wedge_outer_radius, this_wedge_inner_sphericity,
      this_wedge_outer_sphericity, orientation_map_, direction_);
}

double Wedge::operator()(const std::array<double, 3>& source_coords) const {
  return call_impl<double>(source_coords);
}
DataVector Wedge::operator()(
    const std::array<DataVector, 3>& source_coords) const {
  return call_impl<DataVector>(source_coords);
}

std::optional<double> Wedge::original_radius_over_radius(
    const std::array<double, 3>& target_coords, double distorted_radius) const {
  const double radius = magnitude(target_coords);
  CAPTURE_FOR_ERROR(target_coords);
  CAPTURE_FOR_ERROR(distorted_radius);
  CAPTURE_FOR_ERROR(radius);
  CAPTURE_FOR_ERROR(orientation_map_);
  CAPTURE_FOR_ERROR(direction_);

  // First we check the extremal case because if the point is inside the
  // innermost radius surface (the length of the inner cube) or outside the
  // outermost radius, then we know for sure the point isn't in this wedge/block
  if (radius + eps_ < this_wedge_inner_surface_.radius / sqrt(3.0) or
      radius - eps_ > this_wedge_outer_surface_.radius) {
    // Parallel::printf(
    //     "Radius %.16f outside inner radius %.16f and outer radius %.16f\n",
    //     radius, inner_surface_.radius, outer_surface_.radius);
    return std::nullopt;
  }

  // Couple protections that would make a point completely outside of the domain
  // of validity for any wedge
  if (equal_within_roundoff(radius, 0.0) or
      equal_within_roundoff(distorted_radius, 1.0) or distorted_radius > 1.0) {
    Parallel::printf(
        "Wedge transition: Can't get original radius. Target radius = %.16f, "
        "distorted radius = %.16f\n",
        radius, distorted_radius);
    return std::nullopt;
  }

  // Before calculating inner/outer distance, we much check if the angle of the
  // point is actually in this block. We know the opening angle of the wedge is
  // pi/2, so we check the angles in the x and y direction (by removing the
  // other coordinate from the radius).
  const double z_squared = square(target_coords[2]);
  const double radius_xz = sqrt(square(target_coords[0]) + z_squared);
  const double radius_yz = sqrt(square(target_coords[1]) + z_squared);
  // To avoid dividing by zero, if the point has x=0,z=0, then the angle in the
  // x direction is 0 so cos = 1. Same for y direction.
  const double cos_theta_x =
      radius_xz == 0.0 ? 1.0 : target_coords[2] / radius_xz;
  const double cos_theta_y =
      radius_yz == 0.0 ? 1.0 : target_coords[2] / radius_yz;
  const double root_2_over_2 = 0.5 * M_SQRT2;

  CAPTURE_FOR_ERROR(radius_xz);
  CAPTURE_FOR_ERROR(radius_yz);
  CAPTURE_FOR_ERROR(cos_theta_x);
  CAPTURE_FOR_ERROR(cos_theta_y);

  // Parallel::printf(
  //     "In Wedge transition:\n"
  //     " Point: %s\n"
  //     " Radius: %.16f\n"
  //     " Distorted radius: %.16f\n"
  //     " Direction: %s\n",
  //     target_coords, radius, distorted_radius, direction_);

  if (cos_theta_x >= root_2_over_2) {
    // Point is in one of +z,+y,-y wedge. Check cos_theta_y now.
    if (cos_theta_y >= root_2_over_2) {
      // Point is in +z wedge
      if (direction_ != Direction<3>::upper_zeta()) {
        // Parallel::printf(" Not +z, 0\n");
        return std::nullopt;
      }
    } else if (target_coords[1] > 0.0) {
      // Point is in +y wedge
      if (direction_ != Direction<3>::upper_eta()) {
        // Parallel::printf(" Not +y, 1\n");
        return std::nullopt;
      }
    } else if (target_coords[1] < 0.0) {
      // Point is in -y wedge
      if (direction_ != Direction<3>::lower_eta()) {
        // Parallel::printf(" Not -y, 2\n");
        return std::nullopt;
      }
    } else {
      ERROR("Inconsistent angle +z,+y,-y");
    }
  } else if (cos_theta_x <= -root_2_over_2) {
    // Point is in one of -z,+y,-y wedge. Check cos_theta_y now.
    if (cos_theta_y <= -root_2_over_2) {
      // Point is in -z wedge
      if (direction_ != Direction<3>::lower_zeta()) {
        // Parallel::printf(" Not -z, 3\n");
        return std::nullopt;
      }
    } else if (target_coords[1] > 0.0) {
      // Point is in +y wedge
      if (direction_ != Direction<3>::upper_eta()) {
        // Parallel::printf(" Not +y, 4\n");
        return std::nullopt;
      }
    } else if (target_coords[1] < 0.0) {
      // Point is in -y wedge
      if (direction_ != Direction<3>::lower_eta()) {
        // Parallel::printf(" Not -y, 5\n");
        return std::nullopt;
      }
    } else {
      ERROR("Inconsistent angle -z,+y,-y");
    }
  } else if (cos_theta_x < root_2_over_2 and cos_theta_x > -root_2_over_2) {
    // Point is in one of +y,-y,+x,-x. Check phi_z now.
    const double phi_z = atan2(target_coords[1], target_coords[0]);
    const double three_pi_over_four = 3.0 * M_PI_4;
    if (phi_z >= -M_PI_4 and phi_z <= M_PI_4) {
      // Point is in +x wedge
      if (direction_ != Direction<3>::upper_xi()) {
        // Parallel::printf(" Not +x, 6\n");
        return std::nullopt;
      }
    } else if (phi_z >= three_pi_over_four or phi_z <= -three_pi_over_four) {
      // Point is in -x wedge
      if (direction_ != Direction<3>::lower_xi()) {
        // Parallel::printf(" Not -x, 7\n");
        return std::nullopt;
      }
    } else if (target_coords[1] > 0.0) {
      // Point is in +y wedge
      if (direction_ != Direction<3>::upper_eta()) {
        // Parallel::printf(" Not +y, 8\n");
        return std::nullopt;
      }
    } else if (target_coords[1] < 0.0) {
      // Point is in -y wedge
      if (direction_ != Direction<3>::lower_eta()) {
        // Parallel::printf(" Not -y, 9\n");
        return std::nullopt;
      }
    } else {
      ERROR("Inconsistent angle +y,-y,+x,-x");
    }
  } else {
    ERROR("Inconsistent angles.");
  }

  // We are in the proper wedge direction. It is safe to assume the coordinates
  // can be rotated from their original orientation to now "aligned" with +z.
  // Therefore we can continue with computing the original radius.

  const std::array<double, 3> rotated_coords =
      discrete_rotation(orientation_map_.inverse_map(), target_coords);
  CAPTURE_FOR_ERROR(rotated_coords);
  // We use the overall surfaces here because we care about the distances with
  // respect to the innermost and outermost surface as that is what the linear
  // falloff is between.
  const double outer_distance = overall_outer_surface_.distance(rotated_coords);
  const double inner_distance = overall_inner_surface_.distance(rotated_coords);
  const double distance_difference = outer_distance - inner_distance;

  // Now we need to check if the point is actually in this wedge and not a
  // different concentric wedge. We do this by checking the distances
  if (radius + eps_ < inner_distance or radius - eps_ > outer_distance) {
    Parallel::printf(
        "Distance %.16f outside inner distance %.16f and outer distance "
        "%.16f\n",
        radius, inner_distance, outer_distance);
    return std::nullopt;
  }

  // Now we can be sure that this point is in this wedge.

  // If distorted radius is 0, this means the map is the identity so the radius
  // and the original radius are equal. Also we don't want to divide by 0 below
  if (equal_within_roundoff(distorted_radius, 0.0)) {
    return std::optional{1.0};
  }

  // If we are at the outer distance, then the transition function is 0 so the
  // map is again the identity so the radius and original radius are equal
  if (equal_within_roundoff(radius, outer_distance)) {
    return std::optional{1.0};
  }

  // If we are at the inner distance, then the transition function is 1 so we
  // don't need to do a quadratic find. We have r/rtil = 1/(1-SumYlm)
  if (equal_within_roundoff(radius, inner_distance)) {
    return std::optional{1.0 / (1.0 - distorted_radius)};
  }

  // TODO: Comment functional form
  const double a = radius;
  const double c = -distance_difference / distorted_radius;
  const double b = -c - outer_distance;

  // QUESTION: Is this the right one to use?
  std::optional<std::array<double, 2>> roots = real_roots(a, b, c);

  if (roots.has_value()) {
    for (const double root : roots.value()) {
      // Max value is radius of outer sphere (regardless of sphericity) divided
      // by length of inner cube (again regardless of sphericity). Any value
      // bigger than this is not possible
      if (root > 0.0 and
          root < this_wedge_outer_surface_.radius /
                     (this_wedge_inner_surface_.radius / sqrt(3.0))) {
        return std::optional{root};
      }
    }
    using ::operator<<;
    Parallel::printf(
        "No roots in range:\n"
        " a = %.16f\n"
        " b = %.16f\n"
        " c = %.16f\n"
        " Roots: %s\n",
        a, b, c, roots);
  } else {
    using ::operator<<;
    Parallel::printf(
        "Could not compute r/rtil:\n"
        " Target coords = %s\n"
        " Target radius = %.16f\n"
        " Inner distance = %.16f\n"
        " Outer distance = %.16f\n"
        " a = %.16f\n"
        " b = %.16f\n"
        " c = %.16f\n",
        target_coords, radius, inner_distance, outer_distance, a, b, c);
    return std::nullopt;
  }
}

double Wedge::map_over_radius(
    const std::array<double, 3>& source_coords) const {
  return map_over_radius_impl<double>(source_coords);
}
DataVector Wedge::map_over_radius(
    const std::array<DataVector, 3>& source_coords) const {
  return map_over_radius_impl<DataVector>(source_coords);
}

std::array<double, 3> Wedge::gradient(
    const std::array<double, 3>& source_coords) const {
  return gradient_impl<double>(source_coords);
}
std::array<DataVector, 3> Wedge::gradient(
    const std::array<DataVector, 3>& source_coords) const {
  return gradient_impl<DataVector>(source_coords);
}

template <typename T>
T Wedge::call_impl(const std::array<T, 3>& source_coords) const {
  const std::array<T, 3> rotated_coords =
      discrete_rotation(orientation_map_.inverse_map(), source_coords);
  check_distances(rotated_coords);
  T outer_distance = overall_outer_surface_.distance(rotated_coords);

  return (outer_distance - magnitude(rotated_coords)) /
         (outer_distance - overall_inner_surface_.distance(rotated_coords));
}

template <typename T>
T Wedge::map_over_radius_impl(const std::array<T, 3>& source_coords) const {
  const std::array<T, 3> rotated_coords =
      discrete_rotation(orientation_map_.inverse_map(), source_coords);
  check_distances(rotated_coords);
  const T radius = magnitude(rotated_coords);
  const T outer_distance = overall_outer_surface_.distance(rotated_coords);

  return (outer_distance - radius) /
         ((outer_distance - overall_inner_surface_.distance(rotated_coords)) *
          radius);
}

template <typename T>
std::array<T, 3> Wedge::gradient_impl(
    const std::array<T, 3>& source_coords) const {
  // If both surfaces are spherical then we short circuit because the distances
  // are constant and we only need to take a derivative of r.
  // (grad f)_i = -(x_i/r)/(D_out - D_in)
  const std::array<T, 3> rotated_coords =
      discrete_rotation(orientation_map_.inverse_map(), source_coords);
  check_distances(rotated_coords);
  if (overall_inner_surface_.sphericity == 1.0 and
      overall_outer_surface_.sphericity == 1.0) {
    const T one_over_denom =
        1.0 / (magnitude(rotated_coords) *
               (overall_outer_surface_.distance(rotated_coords) -
                overall_inner_surface_.distance(rotated_coords)));

    return discrete_rotation(orientation_map_,
                             -rotated_coords * one_over_denom);
  }

  const T radius = magnitude(rotated_coords);

  const T one_over_radius = 1.0 / radius;
  T outer_distance = overall_outer_surface_.distance(rotated_coords);
  const T one_over_denom =
      1.0 / (outer_distance - overall_inner_surface_.distance(rotated_coords));
  T& outer_distance_minus_radius = outer_distance;
  outer_distance_minus_radius = outer_distance - radius;

  // Avoid roundoff if we are at outer boundary
  for (size_t i = 0; i < get_size(radius); i++) {
    if (equal_within_roundoff(get_element(outer_distance_minus_radius, i),
                              0.0)) {
      get_element(outer_distance_minus_radius, i) = 0.0;
    }
  }

  // Regardless of the sphericities below, we always need this factor in the
  // first term so we calculate it now.
  std::array<T, 3> result = -1.0 * rotated_coords * one_over_radius;

  const auto make_factor = [&one_over_radius](const Surface& surface) -> T {
    return (1.0 - surface.sphericity) * surface.radius * cube(one_over_radius) /
           sqrt(3.0);
  };

  // We can make some simplifications if either of the surfaces are spherical
  // because then the derivative of the distance is zero since it's constant. In
  // the first two branches, it's safe to assume the other sphericity isn't 1
  // because of the above check.
  // TODO: Add more comments below
  if (overall_outer_surface_.sphericity == 1.0) {
    const T inner_surface_factor = make_factor(overall_inner_surface_);

    for (size_t i = 0; i < 2; i++) {
      gsl::at(result, i) -= outer_distance_minus_radius * inner_surface_factor *
                            rotated_coords[2] * gsl::at(rotated_coords, i) *
                            one_over_denom;
    }

    result[2] += outer_distance_minus_radius * inner_surface_factor *
                 (square(radius) - square(rotated_coords[2])) * one_over_denom;
  } else if (overall_inner_surface_.sphericity == 1.0) {
    std::array<T, 3> outer_deriv_distance =
        make_array<3>(make_factor(overall_outer_surface_));
    for (size_t i = 0; i < 2; i++) {
      gsl::at(outer_deriv_distance, i) *=
          -1.0 * rotated_coords[2] * gsl::at(rotated_coords, i);
    }
    outer_deriv_distance[2] *= (square(radius) - square(rotated_coords[2]));

    result += outer_deriv_distance *
              (1.0 - outer_distance_minus_radius * one_over_denom);
  } else {
    const T inner_surface_factor = make_factor(overall_inner_surface_);
    std::array<T, 3> outer_deriv_distance =
        make_array<3>(make_factor(overall_outer_surface_));
    for (size_t i = 0; i < 2; i++) {
      gsl::at(outer_deriv_distance, i) *=
          -1.0 * rotated_coords[2] * gsl::at(rotated_coords, i);
    }
    outer_deriv_distance[2] *= (square(radius) - square(rotated_coords[2]));

    result += outer_deriv_distance;
    for (size_t i = 0; i < 3; i++) {
      if (i != 2) {
        gsl::at(result, i) -= outer_distance_minus_radius *
                              (gsl::at(outer_deriv_distance, i) -
                               (inner_surface_factor * rotated_coords[2] *
                                gsl::at(rotated_coords, i))) *
                              one_over_denom;
      } else {
        gsl::at(result, i) -= outer_distance_minus_radius *
                              (gsl::at(outer_deriv_distance, i) -
                               (inner_surface_factor *
                                (square(radius) - square(rotated_coords[2])))) *
                              one_over_denom;
      }
    }
  }

  // Finally, need one more factor of D_out - D_in in the denominator and to
  // rotate it back to the proper orientation
  return discrete_rotation(orientation_map_, result * one_over_denom);
}

bool Wedge::operator==(const ShapeMapTransitionFunction& other) const {
  if (dynamic_cast<const Wedge*>(&other) == nullptr) {
    return false;
  }
  const Wedge& other_ref = *dynamic_cast<const Wedge*>(&other);
  return overall_inner_surface_ == other_ref.overall_inner_surface_ and
         overall_outer_surface_ == other_ref.overall_outer_surface_ and
         this_wedge_inner_surface_ == other_ref.this_wedge_inner_surface_ and
         this_wedge_outer_surface_ == other_ref.this_wedge_outer_surface_ and
         orientation_map_ == other_ref.orientation_map_;
}

bool Wedge::operator!=(const ShapeMapTransitionFunction& other) const {
  return not(*this == other);
}

// checks that the magnitudes are all between `r_min_` and `r_max_`
template <typename T>
void Wedge::check_distances([
    [maybe_unused]] const std::array<T, 3>& coords) const {
#ifdef SPECTRE_DEBUG
  const T mag = magnitude(coords);
  const T inner_distance = inner_surface_.distance(coords);
  const T outer_distance = outer_surface_.distance(coords);
  for (size_t i = 0; i < get_size(mag); ++i) {
    if (get_element(mag, i) + eps_ < get_element(inner_distance, i) or
        get_element(mag, i) - eps_ > get_element(outer_distance, i)) {
      ERROR(
          "The Wedge transition map was called with coordinates outside "
          "the set inner and outer surfaces. The inner radius and sphericity "
          "are (r="
          << inner_surface_.radius << ",s=" << inner_surface_.sphericity
          << ") and the outer radius and sphericity are (r="
          << outer_surface_.radius << ",s=" << outer_surface_.sphericity
          << "). The inner distance is " << get_element(inner_distance, i)
          << ", the outer distance is " << get_element(outer_distance, i)
          << ". The requested point has radius: " << get_element(mag, i));
    }
  }
#endif  // SPECTRE_DEBUG
}

void Wedge::pup(PUP::er& p) {
  ShapeMapTransitionFunction::pup(p);
  size_t version = 0;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | overall_inner_surface_;
    p | overall_outer_surface_;
    p | this_wedge_inner_surface_;
    p | this_wedge_outer_surface_;
    p | orientation_map_;
    p | direction_;
  }
}

Wedge::Wedge(CkMigrateMessage* const msg) : ShapeMapTransitionFunction(msg) {}

PUP::able::PUP_ID Wedge::my_PUP_ID = 0;

}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
