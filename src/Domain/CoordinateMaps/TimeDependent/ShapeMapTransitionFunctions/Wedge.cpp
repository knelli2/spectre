// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Wedge.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <pup.h>
#include <pup_stl.h>

#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Falloff.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/Structure/Direction.hpp"
#include "NumericalAlgorithms/RootFinding/QuadraticEquation.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/CaptureForError.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {
template <typename T>
T Wedge::Surface::distance(const std::array<T, 3>& coords,
                           const std::optional<size_t>& axis) const {
  // Short circuit if it's a sphere. Then the distance is trivially the radius
  // of this surface
  if (sphericity == 1.0) {
    return make_with_value<T>(coords[0], radius);
  }

  T non_spherical_contribution{};

  if (axis.has_value()) {
#ifdef SPECTRE_DEBUG
    for (size_t i = 0; i < get_size(coords[0]); i++) {
      if (get_element(gsl::at(coords, axis.value()), i) == 0.0) {
        ERROR(
            "The Wedge transition map was called with coordinates outside of "
            "its wedge. The "
            << (axis.value() == 0_st ? "x" : (axis.value() == 1_st ? "y" : "z"))
            << "-coordinate shouldn't be 0 but it is.");
      }
    }
#endif
    non_spherical_contribution = (1.0 - sphericity) / sqrt(3.0) *
                                 magnitude(coords) /
                                 abs(gsl::at(coords, axis.value()));
  } else {
    non_spherical_contribution =
        (1.0 - sphericity) / sqrt(3.0) * magnitude(coords) /
        blaze::max(abs(coords[0]), abs(coords[1]), abs(coords[2]));
  }

  return radius * (non_spherical_contribution + sphericity);
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

Wedge::Wedge(const double inner_radius, const double outer_radius,
             const double inner_sphericity, const double outer_sphericity,
             const Falloff falloff, const size_t axis)
    : inner_surface_(Surface{inner_radius, inner_sphericity}),
      outer_surface_(Surface{outer_radius, outer_sphericity}),
      falloff_(falloff),
      axis_(axis) {
  if (axis_ > 2) {
    ERROR_NO_TRACE(
        "The wedge transition function can only be constructed with an axis of "
        "0, 1, or 2. Not "
        << axis_);
  }
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

  // Couple protections that would make a point completely outside of the domain
  // of validity for any wedge
  if (equal_within_roundoff(radius, 0.0) or
      equal_within_roundoff(distorted_radius, 1.0) or distorted_radius > 1.0) {
    return std::nullopt;
  }

  // First we check the extremal case of being outside the outermost radius. We
  // can check the outermost radius because its surface doesn't move. We can't
  // check the innermost surface because the inner bound is the origin which we
  // already checked above.
  if (radius - eps_ > outer_surface_.radius) {
    return std::nullopt;
  }

  // If distorted radius is 0, this means the map is the identity so the radius
  // and the original radius are equal. Also we don't want to divide by 0 below
  if (equal_within_roundoff(distorted_radius, 0.0)) {
    return std::optional{1.0};
  }

  const double outer_distance = outer_surface_.distance(target_coords);
  const double inner_distance = inner_surface_.distance(target_coords);
  const double distance_difference = outer_distance - inner_distance;

  // If we are at the overall outer distance, then the transition function is 0
  // so the map is again the identity so the radius and original radius are
  // equal. We can't check the overall inner distance because that has been
  // distorted so we don't know where the mapped inner distance is.
  if (equal_within_roundoff(radius, outer_distance)) {
    return std::optional{1.0};
  }

  if (falloff_ == Falloff::Linear) {
    // Solving equation of the form rtil*x^2 + ((D_o - D_i)/SumYlm - D_o)*x -
    // (D_o
    // - D_i)/SumYlm = 0
    const double a = radius;
    const double c = -distance_difference / distorted_radius;
    const double b = -c - outer_distance;

    std::optional<std::array<double, 2>> roots = real_roots(a, b, c);

    if (roots.has_value()) {
      for (const double root : roots.value()) {
        // Check if the root is positive and within the inner and outer distance
        // (divided by the mapped radius)
        if (root > 0.0 and root + eps_ >= inner_distance / radius and
            root - eps_ <= outer_distance / radius) {
          return std::optional{root};
        }
      }
      return std::nullopt;
    } else {
      return std::nullopt;
    }
  } else {
    const double a = -inner_distance / distance_difference;
    const double b = -a * outer_distance;

    const double denom = 1.0 - distorted_radius * a;

    if (equal_within_roundoff(denom, 0.0)) {
      return std::nullopt;
    }

    const double original_radius = (radius + distorted_radius * b) / denom;

    if (original_radius + eps_ >= inner_distance and
        original_radius - eps_ <= outer_distance) {
      return std::optional{original_radius / radius};
    } else {
      return std::nullopt;
    }
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
  check_distances(source_coords);
  const T radius = magnitude(source_coords);
  const T outer_distance = outer_surface_.distance(source_coords, axis_);
  const T inner_distance = inner_surface_.distance(source_coords, axis_);
  const T one_over_denom = 1.0 / (outer_distance - inner_distance);

  if (falloff_ == Falloff::Linear) {
    return (outer_distance - radius) * one_over_denom;
  } else {
    const T a = -inner_distance * one_over_denom;
    const T b = -a * outer_distance;

    return a + b / radius;
  }
}

template <typename T>
T Wedge::map_over_radius_impl(const std::array<T, 3>& source_coords) const {
  check_distances(source_coords);
  const T radius = magnitude(source_coords);
  const T outer_distance = outer_surface_.distance(source_coords, axis_);
  const T inner_distance = inner_surface_.distance(source_coords, axis_);
  const T one_over_denom = 1.0 / (outer_distance - inner_distance);

  if (falloff_ == Falloff::Linear) {
    return (outer_distance - radius) * one_over_denom / radius;
  } else {
    const T a_over_r = -inner_distance * one_over_denom / radius;
    const T b_over_r = -a_over_r * outer_distance;

    return a_over_r + b_over_r / radius;
  }
}

template <typename T>
std::array<T, 3> Wedge::gradient_impl(
    const std::array<T, 3>& source_coords) const {
  // If both surfaces are spherical then we short circuit because the distances
  // are constant and we only need to take a derivative of r.
  // (grad f)_i = -(x_i/r)/(D_out - D_in)
  check_distances(source_coords);

  const T radius = magnitude(source_coords);
  const T one_over_radius = 1.0 / radius;

  const T outer_distance = outer_surface_.distance(source_coords, axis_);
  const T inner_distance = inner_surface_.distance(source_coords, axis_);
  const T one_over_distance_difference =
      1.0 / (outer_distance - inner_distance);

  CAPTURE_FOR_ERROR(source_coords);
  CAPTURE_FOR_ERROR(radius);
  CAPTURE_FOR_ERROR(outer_distance);
  CAPTURE_FOR_ERROR(inner_distance);

  // Special case for both spherical boundaries
  if (inner_surface_.sphericity == 1.0 and outer_surface_.sphericity == 1.0) {
    if (falloff_ == Falloff::Linear) {
      return -1.0 * source_coords * one_over_radius *
             one_over_distance_difference;
    } else {
      const T b =
          outer_distance * inner_distance * one_over_distance_difference;
      return -1.0 * b * source_coords * cube(one_over_radius);
    }
  }

  const auto surface_gradient = [&](const Surface& surface) {
    if (surface.sphericity == 1.0) {
      return make_array<3, T>(make_with_value<T>(source_coords[0], 0.0));
    }

    const size_t axis_plus_one = (axis_ + 1) % 3;
    const size_t axis_plus_two = (axis_ + 2) % 3;

    const T& axis_coord = gsl::at(source_coords, axis_);

    const double factor =
        surface.radius * (1.0 - surface.sphericity) / sqrt(3.0);

    std::array<T, 3> grad =
        make_array<3, T>(factor * one_over_radius / abs(axis_coord));

    // Dividing by axis_coord here takes care of the sgn(axis_coord) that would
    // have been necessary
    gsl::at(grad, axis_) *= -(square(radius) - square(axis_coord)) / axis_coord;
    gsl::at(grad, axis_plus_one) *= gsl::at(source_coords, axis_plus_one);
    gsl::at(grad, axis_plus_two) *= gsl::at(source_coords, axis_plus_two);

    return grad;
  };

  const std::array<T, 3> outer_gradient = surface_gradient(outer_surface_);
  const std::array<T, 3> inner_gradient = surface_gradient(inner_surface_);

  return outer_gradient;

  if (falloff_ == Falloff::Linear) {
    return (outer_gradient - source_coords * one_over_radius -
            (outer_distance - radius) * (outer_gradient - inner_gradient) *
                one_over_distance_difference) *
           one_over_distance_difference;
  } else {
    const T a = -inner_distance * one_over_distance_difference;
    const T b = -a * outer_distance;

    const std::array<T, 3> a_gradient =
        (inner_distance * outer_gradient - outer_distance * inner_gradient) *
        square(one_over_distance_difference);

    return a_gradient -
           one_over_radius *
               (outer_gradient * a + outer_distance * a_gradient) -
           b * source_coords * cube(one_over_radius);
  }
}

bool Wedge::operator==(const ShapeMapTransitionFunction& other) const {
  if (dynamic_cast<const Wedge*>(&other) == nullptr) {
    return false;
  }
  const Wedge& other_ref = *dynamic_cast<const Wedge*>(&other);
  return inner_surface_ == other_ref.inner_surface_ and
         outer_surface_ == other_ref.outer_surface_ and
         falloff_ == other_ref.falloff_;
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
  const T inner_distance = inner_surface_.distance(coords, axis_);
  const T outer_distance = outer_surface_.distance(coords, axis_);
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
    p | inner_surface_;
    p | outer_surface_;
    p | falloff_;
    p | axis_;
  }
}

Wedge::Wedge(CkMigrateMessage* const msg) : ShapeMapTransitionFunction(msg) {}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
PUP::able::PUP_ID Wedge::my_PUP_ID = 0;

}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
