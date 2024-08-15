// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Wedge.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <pup.h>
#include <pup_stl.h>
#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/Structure/Direction.hpp"
#include "NumericalAlgorithms/RootFinding/QuadraticEquation.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/CaptureForError.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/StdHelpers.hpp"

namespace {
template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));  // NOLINT
}
}  // namespace

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {
Wedge::Surface::Surface(const std::array<double, 3>& center_in,
                        const double radius_in, const double sphericity_in)
    : center(center_in),
      radius(radius_in),
      sphericity(sphericity_in),
      half_cube_length(radius / sqrt(3.0)) {}

void Wedge::Surface::pup(PUP::er& p) {
  p | center;
  p | radius;
  p | sphericity;
  p | half_cube_length;
}

size_t Wedge::axis_index(const Axis axis) {
  return static_cast<size_t>(abs(static_cast<int>(axis)) - 1);
}

double Wedge::axis_sgn(const Axis axis) {
  return static_cast<double>(sgn(static_cast<int>(axis)));
}

std::ostream& operator<<(std::ostream& os, const Wedge::Axis axis) {
  os << (Wedge::axis_sgn(axis) < 0.0 ? "-" : "+") << Wedge::axis_index(axis);
  return os;
}

Wedge::Wedge(const std::array<double, 3>& inner_center,
             const double inner_radius, const double inner_sphericity,
             const std::array<double, 3>& outer_center,
             const double outer_radius, const double outer_sphericity,
             const Axis axis)
    : inner_surface_(inner_center, inner_radius, inner_sphericity),
      outer_surface_(outer_center, outer_radius, outer_sphericity),
      projection_center_(inner_surface_.center - outer_surface_.center),
      axis_(axis) {
  if (projection_center_ != std::array{0.0, 0.0, 0.0} and
      inner_surface_.sphericity != 1.0) {
    ERROR(
        "The sphericity of the inner surface must be exactly 1.0 if the "
        "centers of the inner and outer object are different.\ninner center: "
        << inner_surface_.center
        << "\nouter center: " << outer_surface_.center);
  }
  if (axis_ == Axis::None) {
    ERROR(
        "Axis for Wedge shape map transition function cannot be 'None'. Please "
        "choose another.");
  }
}

template <typename T>
std::array<T, 3> Wedge::compute_inner_surface_vector(
    const std::array<T, 3>& centered_coords, const T& centered_coords_magnitude,
    const std::optional<Axis>& potential_axis) const {
  // If the surface is a sphere, then the vector is trivially the radial vector
  if (inner_surface_.sphericity == 1.0) {
    return inner_surface_.radius * centered_coords / centered_coords_magnitude;
  }

  T non_spherical_contribution =
      (1.0 - inner_surface_.sphericity) / sqrt(3.0) * centered_coords_magnitude;
  // Do some debug checks so we don't divide by zero below if we are given an
  // axis
  if (potential_axis.has_value()) {
    const int index = axis_index(potential_axis.value());
#ifdef SPECTRE_DEBUG
    for (size_t i = 0; i < get_size(centered_coords[0]); i++) {
      if (get_element(gsl::at(centered_coords, index), i) == 0.0) {
        ERROR(
            "The Wedge transition map was called with coordinates outside of "
            "its wedge. The "
            << index << "-coordinate shouldn't be 0 but it is. Coordinate: ("
            << get_element(centered_coords[0], i) << ","
            << get_element(centered_coords[1], i) << ","
            << get_element(centered_coords[2], i) << ")");
      }
    }
#endif
    non_spherical_contribution /= abs(gsl::at(centered_coords, index));
  } else {
    non_spherical_contribution /=
        blaze::max(abs(centered_coords[0]), abs(centered_coords[1]),
                   abs(centered_coords[2]));
  }

  return inner_surface_.radius *
         (non_spherical_contribution + inner_surface_.sphericity) *
         centered_coords / centered_coords_magnitude;
}

template <typename T>
std::array<T, 3> Wedge::compute_outer_surface_vector(
    const std::array<T, 3>& centered_coords, const T& lambda) const {
  return lambda * centered_coords;
}

template <typename T>
T Wedge::lambda_cube(const std::array<T, 3>& centered_coords,
                     const std::optional<Axis>& potential_axis) const {
  const auto lambda_axis = [&centered_coords, this](const Axis axis) -> T {
    const size_t axis_idx = axis_index(axis);
    const double factor = (axis_sgn(axis) * outer_surface_.half_cube_length -
                           gsl::at(projection_center_, axis_idx));
    // Avoid dividing by zero. If a coordinate is zero, then that point can't
    // exist in the wedge whose coordinate is zero. Since we take the min of all
    // positive lambdas to find the correct one, we overwrite this lambda with
    // -max so it'll never be chosen
    T result = make_with_value<T>(centered_coords[0], 0.0);
    for (size_t i = 0; i < get_size(result); i++) {
      if (get_element(gsl::at(centered_coords, axis_idx), i) == 0.0) {
        get_element(result, i) = -std::numeric_limits<double>::max();
      } else {
        get_element(result, i) =
            factor / get_element(gsl::at(centered_coords, axis_idx), i);
      }
    }
    // Parallel::printf(
    //     "lambda cube:\n"
    //     " axis: %s\n"
    //     " lambda: %s\n",
    //     axis, MakeString{} << result);
    return result;
  };

  // If we specified an axis that we are on, then use that lambda. This saves
  // some computation.
  if (potential_axis.has_value()) {
    return lambda_axis(potential_axis.value());
  }

  // Otherwise we need to loop over all of them
  const std::array all_axes{Axis::PlusZ,  Axis::MinusZ, Axis::PlusY,
                            Axis::MinusY, Axis::PlusX,  Axis::MinusX};
  std::array<T, 6> candidate_lambdas{};

  for (size_t i = 0; i < all_axes.size(); i++) {
    gsl::at(candidate_lambdas, i) = lambda_axis(gsl::at(all_axes, i));
  }

  // We need to take the min of all the positive lambdas
  const double fill_value = std::numeric_limits<double>::max();
  T result = make_with_value<T>(candidate_lambdas[0], fill_value);
  std::vector<Axis> chosen_axes(get_size(result), Axis::None);
  for (size_t i = 0; i < get_size(result); i++) {
    for (size_t j = 0; j < candidate_lambdas.size(); j++) {
      if (get_element(gsl::at(candidate_lambdas, j), i) > 0.0 and
          get_element(gsl::at(candidate_lambdas, j), i) <
              get_element(result, i)) {
        get_element(result, i) = get_element(gsl::at(candidate_lambdas, j), i);
        chosen_axes[i] = gsl::at(all_axes, j);
      }
    }
  }

  // Parallel::printf(
  //     "lambda cube:\n"
  //     " chosen axis: %s\n",
  //     MakeString{} << chosen_axes);

  ASSERT(not alg::any_of(chosen_axes,
                         [&](const Axis& a) { return a == Axis::None; }),
         "No candidate cube lambdas in the Wedge shape map transition function "
         "are negative. This is an internal error. Please file an issue.");

  return result;
}

template <typename T>
T Wedge::lambda_sphere(const std::array<T, 3>& centered_coords,
                       const T& centered_coords_magnitude) const {
  // quadratic equation is
  // a lambda^2 + b lambda + c = 0
  //
  // For this case
  // a = |x^i-P^i|^2
  // b = 2(x^i-P^i)(P^j-C^j)\delta_{ij}
  // c = |P^i-C^i|^2 - R^2

  const T a = square(centered_coords_magnitude);
  const T b = 2.0 * dot(centered_coords, projection_center_);
  const T c =
      make_with_value<T>(b, dot(projection_center_, projection_center_) -
                                square(outer_surface_.radius));

  // The root that we are looking for is negative (since a positive root would
  // be for the opposite side of the sphere)
  // if constexpr (std::is_same_v<T, double>) {
  //   Parallel::printf(
  //       "lambda sphere:\n"
  //       " a: %.16f\n"
  //       " b: %.16f\n"
  //       " c: %.16f\n"
  //       " roots: %s\n",
  //       a, b, c, real_roots(a, b, c));
  // }
  return smallest_root_greater_than_value_within_roundoff(
      a, b, c, -std::numeric_limits<double>::epsilon() * 100.0);
}

template <typename T>
T Wedge::compute_lambda(const std::array<T, 3>& centered_coords,
                        const T& centered_coords_magnitude,
                        const std::optional<Axis>& potential_axis) const {
  T result{};
  // Avoid extra computation if possible
  if (outer_surface_.sphericity == 1.0) {
    result = lambda_sphere(centered_coords, centered_coords_magnitude);
  } else if (outer_surface_.sphericity == 0.0) {
    result = lambda_cube(centered_coords, potential_axis);
  } else {
    result = (1.0 - outer_surface_.sphericity) *
                 lambda_cube(centered_coords, potential_axis) +
             outer_surface_.sphericity *
                 lambda_sphere(centered_coords, centered_coords_magnitude);
  }

  return result;
}

template <typename T>
std::array<T, 3> Wedge::lambda_cube_gradient(
    const T& lambda_cube, const std::array<T, 3>& centered_coords) const {
  const size_t axis_idx = axis_index(axis_);

  std::array<T, 3> result =
      make_array<3>(make_with_value<T>(centered_coords[0], 0.0));
  gsl::at(result, axis_idx) = -lambda_cube / gsl::at(centered_coords, axis_idx);

  // Parallel::printf(
  //     "lambda cube gradient:\n"
  //     " lambda cube: %s\n"
  //     " Centered coord: %s\n",
  //     MakeString{} << lambda_cube, centered_coords);

  return result;
}

template <typename T>
std::array<T, 3> Wedge::lambda_sphere_gradient(
    const T& lambda_sphere, const std::array<T, 3>& centered_coords,
    const T& centered_coords_magnitude) const {
  return -1.0 * lambda_sphere *
         (lambda_sphere * centered_coords + projection_center_) /
         (lambda_sphere * square(centered_coords_magnitude) +
          dot(centered_coords, projection_center_));
}

template <typename T>
std::array<T, 3> Wedge::compute_lambda_gradient(
    const std::array<T, 3>& centered_coords,
    const T& centered_coords_magnitude) const {
  std::array<T, 3> result{};
  // Avoid extra computation if possible
  if (outer_surface_.sphericity == 1.0) {
    result = lambda_sphere_gradient(
        lambda_sphere(centered_coords, centered_coords_magnitude),
        centered_coords, centered_coords_magnitude);
  } else if (outer_surface_.sphericity == 0.0) {
    result = lambda_cube_gradient(lambda_cube(centered_coords, {axis_}),
                                  centered_coords);
  } else {
    result = (1.0 - outer_surface_.sphericity) *
                 lambda_cube_gradient(lambda_cube(centered_coords, {axis_}),
                                      centered_coords) +
             outer_surface_.sphericity *
                 lambda_sphere_gradient(
                     lambda_sphere(centered_coords, centered_coords_magnitude),
                     centered_coords, centered_coords_magnitude);
  }

  return result;
}

double Wedge::operator()(const std::array<double, 3>& source_coords) const {
  return call_impl<double>(source_coords);
}

DataVector Wedge::operator()(
    const std::array<DataVector, 3>& source_coords) const {
  return call_impl<DataVector>(source_coords);
}

template <typename T>
T Wedge::call_impl(const std::array<T, 3>& source_coords) const {
  // const std::array<T, 3> centered_coords = source_coords -
  // projection_center_;
  const std::array<T, 3> centered_coords = source_coords;
  const T centered_coords_magnitude = magnitude(centered_coords);

  const T lambda =
      compute_lambda(centered_coords, centered_coords_magnitude, {axis_});

  const T inner_distance =
      inner_surface_.sphericity == 1.0
          ? make_with_value<T>(lambda, inner_surface_.radius)
          : magnitude(compute_inner_surface_vector(
                centered_coords, centered_coords_magnitude, {axis_}));

  const std::array<T, 3> outer_surface_vector =
      compute_outer_surface_vector(centered_coords, lambda);
  const T outer_distance = magnitude(outer_surface_vector);

  // if constexpr (std::is_same_v<T, DataVector>) {
  //   Parallel::printf(
  //       "operator()\n"
  //       " Coord: %s\n"
  //       " Centered coord: %s\n"
  //       " Centered coord mag: %s\n"
  //       " Lambda: %s\n"
  //       " Inner distance: %s\n"
  //       " Outer surface vec: %s\n"
  //       " Outer surface vec mag: %s\n"
  //       " result: %s\n",
  //       source_coords, centered_coords, centered_coords_magnitude, lambda,
  //       inner_distance, outer_surface_vector, outer_distance,
  //       (outer_distance - centered_coords_magnitude) /
  //           (outer_distance - inner_surface_.radius));
  // } else {
  //   Parallel::printf(
  //       "operator()\n"
  //       " Coord: %s\n"
  //       " Centered coord: %s\n"
  //       " Centered coord mag: %.16f\n"
  //       " Lambda: %.16f\n"
  //       " Inner distance: %.16f\n"
  //       " Outer surface vec: %s\n"
  //       " Outer surface vec mag: %.16f\n"
  //       " result: %.16f\n",
  //       source_coords, centered_coords, centered_coords_magnitude, lambda,
  //       inner_distance, outer_surface_vector, outer_distance,
  //       (outer_distance - centered_coords_magnitude) /
  //           (outer_distance - inner_surface_.radius));
  // }

#ifdef SPECTRE_DEBUG
  const T result = (outer_distance - centered_coords_magnitude) /
                   (outer_distance - inner_distance);
  for (size_t i = 0; i < get_size(centered_coords_magnitude); ++i) {
    if (get_element(result, i) + eps_ < 0.0 or
        get_element(result, i) - eps_ > 1.0) {
      const std::array point{get_element(centered_coords[0], i),
                             get_element(centered_coords[1], i),
                             get_element(centered_coords[2], i)};
      ERROR(
          "The Wedge transition map was called with coordinates outside "
          "the set inner and outer surfaces.\nThe requested (centered) point "
          "is "
          << point << "\nThe requested (centered) point has radius "
          << get_element(centered_coords_magnitude, i)
          << "\nThe inner surface has center, "
             "radius, and sphericity (c="
          << inner_surface_.center << ",r=" << inner_surface_.radius
          << ",s=" << inner_surface_.sphericity
          << ")\nThe outer surface has center, radius, and sphericity (c="
          << outer_surface_.center << ",r=" << outer_surface_.radius
          << ",s=" << outer_surface_.sphericity
          << ")\nThe distance to the inner surface is "
          << get_element(inner_distance, i)
          << "\nThe distance to the outer surface is "
          << get_element(outer_distance, i));
    }
  }
#endif

  return (outer_distance - centered_coords_magnitude) /
         (outer_distance - inner_distance);
}

std::optional<double> Wedge::original_radius_over_radius(
    const std::array<double, 3>& target_coords, double distorted_radius) const {
  // const std::array<double, 3> centered_coords =
  //     target_coords - projection_center_;
  const std::array<double, 3> centered_coords = target_coords;
  const double centered_coords_magnitude = magnitude(centered_coords);
  CAPTURE_FOR_ERROR(target_coords);
  CAPTURE_FOR_ERROR(distorted_radius);
  CAPTURE_FOR_ERROR(centered_coords_magnitude);

  // Couple protections that would make a point completely outside of the domain
  // of validity for any wedge
  if (equal_within_roundoff(centered_coords_magnitude, 0.0) or
      equal_within_roundoff(distorted_radius, 1.0) or distorted_radius > 1.0) {
    return std::nullopt;
  }

  const double inner_distance =
      inner_surface_.sphericity == 1.0
          ? inner_surface_.radius
          : magnitude(compute_inner_surface_vector(
                centered_coords, centered_coords_magnitude, std::nullopt));

  const double lambda =
      compute_lambda(centered_coords, centered_coords_magnitude, std::nullopt);

  const double outer_distance =
      magnitude(compute_outer_surface_vector(centered_coords, lambda));

  // First we check the extremal case of being outside the outer distance. We
  // can check the outermost distance because its surface doesn't move. We can't
  // check the innermost surface because the inner bound is the origin which we
  // already checked above.
  if (centered_coords_magnitude - eps_ > outer_distance) {
    return std::nullopt;
  }

  // If distorted radius is 0, this means the map is the identity so the radius
  // and the original radius are equal. Also we don't want to divide by 0 below.
  // We do this check after we check if the point is beyond the outer distance
  // because if a point is outside the distorted frame, this function should
  // return nullopt.
  if (equal_within_roundoff(distorted_radius, 0.0)) {
    return std::optional{1.0};
  }

  // If we are at the overall outer distance, then the transition function is 0
  // so the map is again the identity so the radius and original radius are
  // equal. We can't check the overall inner distance because that has been
  // distorted so we don't know where the mapped inner distance is.
  if (equal_within_roundoff(centered_coords_magnitude, outer_distance)) {
    return std::optional{1.0};
  }

  const double denom = outer_distance - inner_distance + distorted_radius;
  // prevent zero division
  if (equal_within_roundoff(denom, 0.)) {
    return std::nullopt;
  }

  const double original_radius_over_radius =
      (outer_distance * (1.0 + distorted_radius / centered_coords_magnitude) -
       inner_distance) /
      denom;
  const double original_radius =
      original_radius_over_radius * centered_coords_magnitude;

  // Parallel::printf(
  //     "original radius over radius:\n"
  //     " axis: %s\n"
  //     " target coords: %s\n"
  //     " centered target coords: %s\n"
  //     " centered target coords mag: %.16f\n"
  //     " distorted radius: %.16f\n"
  //     " inner distance: %.16f\n"
  //     " lambda: %.16f\n"
  //     " outer distance: %.16f\n"
  //     " original radius: %.16f\n"
  //     " original radius over radius: %.16f\n"
  //     " denom: %.16f\n",
  //     MakeString{} << axis_, target_coords, centered_coords,
  //     centered_coords_magnitude, distorted_radius, inner_distance, lambda,
  //     outer_distance, original_radius, original_radius_over_radius, denom);

  return (original_radius + eps_) >= inner_distance and
                 (original_radius - eps_) <= outer_distance
             ? std::optional<double>{original_radius_over_radius}
             : std::nullopt;
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
std::array<T, 3> Wedge::gradient_impl(
    const std::array<T, 3>& source_coords) const {
// The operator does a debug check of the source coordinates. The extra
// computation isn't a big deal in debug mode and we don't use the result
#ifdef SPECTRE_DEBUG
  operator()(source_coords);
#endif

  // const std::array<T, 3> centered_coords = source_coords -
  // projection_center_;
  const std::array<T, 3> centered_coords = source_coords;
  const T centered_coords_magnitude = magnitude(centered_coords);
  const T one_over_centered_coords_magnitude = 1.0 / centered_coords_magnitude;

  const T lambda =
      compute_lambda(centered_coords, centered_coords_magnitude, {axis_});

  const T inner_distance =
      inner_surface_.sphericity == 1.0
          ? make_with_value<T>(lambda, inner_surface_.radius)
          : magnitude(compute_inner_surface_vector(
                centered_coords, centered_coords_magnitude, {axis_}));

  const std::array<T, 3> outer_surface_vector =
      compute_outer_surface_vector(centered_coords, lambda);
  const T outer_distance = magnitude(outer_surface_vector);

  // This can only be called if the projection center is 0, otherwise this
  // formula won't work. And if the projection center isn't 0, then we require
  // that the sphericity of the inner surface is 1, so we ASSERT that here as
  // well
  const auto inner_surface_gradient = [&]() -> std::array<T, 3> {
    ASSERT(inner_surface_.sphericity != 1.0,
           "Should not be calculating the inner surface gradient when its "
           "sphericity is 1.0 because the gradient is exactly 0.");

    const size_t axis_idx = axis_index(axis_);
    const size_t axis_plus_one = (axis_idx + 1) % 3;
    const size_t axis_plus_two = (axis_idx + 2) % 3;

    const T& axis_coord = gsl::at(centered_coords, axis_idx);

    const double factor =
        inner_surface_.radius * (1.0 - inner_surface_.sphericity) / sqrt(3.0);

    std::array<T, 3> grad = make_array<3, T>(
        factor * one_over_centered_coords_magnitude / abs(axis_coord));

    // Dividing by axis_coord here takes care of the sgn(axis_coord) that would
    // have been necessary
    gsl::at(grad, axis_idx) *=
        -(square(centered_coords_magnitude) - square(axis_coord)) / axis_coord;
    gsl::at(grad, axis_plus_one) *= gsl::at(centered_coords, axis_plus_one);
    gsl::at(grad, axis_plus_two) *= gsl::at(centered_coords, axis_plus_two);

    return grad;
  };

  // We always need to compute the outer gradient regardless of the projection
  // center or the outer sphericity
  const auto outer_surface_gradient = [&]() -> std::array<T, 3> {
    const std::array<T, 3> lambda_gradient =
        compute_lambda_gradient(centered_coords, centered_coords_magnitude);

    // if constexpr (std::is_same_v<T, DataVector>) {
    //   Parallel::printf(
    //       "outer gradient:\n"
    //       " lambda: %s\n"
    //       " lambda gradient: %s\n"
    //       " centered coords mag: %s\n",
    //       lambda, lambda_gradient, centered_coords_magnitude);
    // } else {
    //   Parallel::printf(
    //       "outer gradient:\n"
    //       " lambda: %.16f\n"
    //       " lambda gradient: %s\n"
    //       " centered coords mag: %.16f\n",
    //       lambda, lambda_gradient, centered_coords_magnitude);
    // }

    return lambda * centered_coords * one_over_centered_coords_magnitude +
           lambda_gradient * centered_coords_magnitude;
  };

  const std::array<T, 3> outer_gradient = outer_surface_gradient();

  const T one_over_distance_difference =
      1.0 / (outer_distance - inner_distance);

  // if constexpr (std::is_same_v<T, DataVector>) {
  //   Parallel::printf(
  //       "gradient:\n"
  //       " Axis: %s\n"
  //       " Coord: %s\n"
  //       " Centered coord: %s\n"
  //       " Centered coord mag: %s\n"
  //       " Lambda: %s\n"
  //       " Inner distance: %s\n"
  //       " Inner gradient: %s\n"
  //       " Outer surface vec: %s\n"
  //       " Outer surface vec mag: %s\n"
  //       " Outer surface grad: %s\n",
  //       axis_, source_coords, centered_coords, centered_coords_magnitude,
  //       lambda, inner_distance,
  //       (inner_surface_.sphericity == 1.0 ? std::array<T, 3>{}
  //                                         : inner_surface_gradient()),
  //       outer_surface_vector, outer_distance, outer_gradient);
  // } else {
  //   Parallel::printf(
  //       "gradient:\n"
  //       " Axis: %s\n"
  //       " Coord: %s\n"
  //       " Centered coord: %s\n"
  //       " Centered coord mag: %.16f\n"
  //       " Lambda: %.16f\n"
  //       " Inner distance: %.16f\n"
  //       " Inner gradient: %s\n"
  //       " Outer surface vec: %s\n"
  //       " Outer surface vec mag: %.16f\n"
  //       " Outer surface grad: %s\n",
  //       axis_, source_coords, centered_coords, centered_coords_magnitude,
  //       lambda, inner_distance,
  //       (inner_surface_.sphericity == 1.0 ? std::array<T, 3>{}
  //                                         : inner_surface_gradient()),
  //       outer_surface_vector, outer_distance, outer_gradient);
  // }

  // Avoid allocating an array of 0 if the inner surface is a sphere
  if (inner_surface_.sphericity == 1.0) {
    return (outer_gradient -
            centered_coords * one_over_centered_coords_magnitude -
            (outer_distance - centered_coords_magnitude) * outer_gradient *
                one_over_distance_difference) *
           one_over_distance_difference;
  } else {
    const std::array<T, 3> inner_gradient = inner_surface_gradient();

    return (outer_gradient -
            centered_coords * one_over_centered_coords_magnitude -
            (outer_distance - centered_coords_magnitude) *
                (outer_gradient - inner_gradient) *
                one_over_distance_difference) *
           one_over_distance_difference;
  }
}

bool Wedge::operator==(const ShapeMapTransitionFunction& other) const {
  if (dynamic_cast<const Wedge*>(&other) == nullptr) {
    return false;
  }
  const auto surface_equal = [](const Surface& s1, const Surface& s2) {
    return s1.center == s2.center and s1.radius == s2.radius and
           s1.sphericity == s2.sphericity and
           s1.half_cube_length == s2.half_cube_length;
  };
  const Wedge& other_ref = *dynamic_cast<const Wedge*>(&other);
  return surface_equal(inner_surface_, other_ref.inner_surface_) and
         surface_equal(outer_surface_, other_ref.outer_surface_) and
         axis_ == other_ref.axis_;
}

bool Wedge::operator!=(const ShapeMapTransitionFunction& other) const {
  return not(*this == other);
}

namespace {
// struct to allow unpacking version 0 of this class
struct LegacySurface {
  double radius{};
  double sphericity{};

  void pup(PUP::er& p) {
    p | radius;
    p | sphericity;
  }
};
}  // namespace

void Wedge::pup(PUP::er& p) {
  ShapeMapTransitionFunction::pup(p);
  size_t version = 1;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 1) {
    p | inner_surface_;
    p | outer_surface_;
    p | projection_center_;
    p | axis_;
  } else if (p.isUnpacking()) {
    LegacySurface inner_surface{};
    LegacySurface outer_surface{};
    size_t axis = 0;
    p | inner_surface;
    p | outer_surface;
    p | axis;

    inner_surface_ = Surface{
        {0.0, 0.0, 0.0},
        inner_surface.radius,
        inner_surface.sphericity,
    };
    outer_surface_ = Surface{
        {0.0, 0.0, 0.0},
        outer_surface.radius,
        outer_surface.sphericity,
    };
    // Unfortunately we can't recover the sign for the axis so we just set it
    // to positive. This is wrong though.
    axis_ = static_cast<Axis>(axis + 1);
    projection_center_ = std::array{0.0, 0.0, 0.0};
  }
}

Wedge::Wedge(CkMigrateMessage* const msg) : ShapeMapTransitionFunction(msg) {}

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
PUP::able::PUP_ID Wedge::my_PUP_ID = 0;
}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
