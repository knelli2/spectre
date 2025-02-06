// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <pup.h>

#include "DataStructures/Blaze/IntegerPow.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/Math.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {

SphereTransition::SphereTransition(const double r_min, const double r_max,
                                   const bool reverse, const bool interior)
    : r_min_(r_min), r_max_(r_max), interior_(interior) {
  if (r_min <= 0.) {
    ERROR("The minimum radius must be greater than 0 but is " << r_min);
  }
  if (r_max <= r_min) {
    ERROR(
        "The maximum radius must be greater than the minimum radius but "
        "r_max =  "
        << r_max << ", and r_min = " << r_min);
  }
  if (interior and reverse) {
    ERROR("Cannot be reverse while also in the interior.");
  }
  a_ = -1.0 / (r_max - r_min);
  b_ = -a_ * r_max;
  if (reverse) {
    a_ *= -1.0;
    b_ = 1.0 - b_;
  }
}

double SphereTransition::operator()(
    const std::array<double, 3>& source_coords,
    const std::optional<size_t>& one_over_radius_power) const {
  const double mag = magnitude(source_coords);

  if (UNLIKELY(one_over_radius_power.has_value() and
               equal_within_roundoff(mag, 0.0))) {
    ERROR("Trying to divide by a point "
          << source_coords
          << " with radius zero in SphereTransition operator.");
  }

  if (interior_) {
    return 1.0 / (r_min_ *
                  integer_pow(mag, static_cast<int>(
                                       one_over_radius_power.value_or(0_st))));
  }

  if (UNLIKELY(mag < r_min_ - eps_)) {
    ERROR("SphereTransition coord " << source_coords << " with radius " << mag
                                    << " is within r_min of " << r_min_
                                    << ", but the class was not constructed "
                                       "for the interior of the sphere.");
  } else if (UNLIKELY(mag > r_max_ + eps_)) {
    ERROR("SphereTransition coord " << source_coords << " with radius " << mag
                                    << "is beyond r_max of " << r_max_ << ".");
  }

  double result = (a_ * mag + b_);
  // Avoid roundoff
  result = std::clamp(result, 0.0, 1.0);

  return result /
         integer_pow(
             mag, static_cast<int>(1 + one_over_radius_power.value_or(0_st)));
}

DataVector SphereTransition::operator()(
    const std::array<DataVector, 3>& source_coords,
    const std::optional<size_t>& one_over_radius_power) const {
  DataVector result = source_coords[0];
  // Go point by point
  for (size_t i = 0; i < source_coords[0].size(); i++) {
    result[i] = (*this)(std::array{source_coords[0][i], source_coords[1][i],
                                   source_coords[2][i]},
                        one_over_radius_power);
  }
  return result;
}

std::optional<double> SphereTransition::original_radius_over_radius(
    const std::array<double, 3>& target_coords,
    double radial_distortion) const {
  const double mag = magnitude(target_coords);
  // If we are at the center, the radius is the same
  if (UNLIKELY(equal_within_roundoff(mag, 0.0))) {
    return interior_ ? std::optional{1.0} : std::nullopt;
  }

  // a_ being positive is a sentinel for reversed.
  // If we aren't reversed, check near or within r_min_.
  if (a_ < 0.0) {
    if (mag + radial_distortion < r_min_ - eps_) {
      return interior_ ? std::optional{r_min_ / (r_min_ - radial_distortion)}
                       : std::nullopt;
    } else if (equal_within_roundoff(mag, r_min_)) {
      return std::optional{r_min_ / (r_min_ - radial_distortion)};
    }
  }

  // Beyond the range of validity for both reversed and not reversed
  if ((a_ < 0.0 and mag > r_max_ + eps_) or
      (a_ > 0.0 and mag < r_min_ - eps_)) {
    return std::nullopt;
  }

  // No distortion means our point is the same
  if (equal_within_roundoff(radial_distortion, 0.0)) {
    return std::optional{1.0};
  }

  // At the f=0 boundary.
  if ((a_ < 0.0 and equal_within_roundoff(mag, r_max_)) or
      (a_ > 0.0 and equal_within_roundoff(mag, r_min_))) {
    return std::optional{1.0};
  }

  const double denom = 1. - radial_distortion * a_;
  // prevent zero division
  if (UNLIKELY(equal_within_roundoff(denom, 0.))) {
    return std::nullopt;
  }

  const double original_radius = (mag + radial_distortion * b_) / denom;

  // Check at or beyond f=1 boundary for reversed
  if (a_ > 0.0) {
    if (equal_within_roundoff(original_radius, r_max_)) {
      return std::optional{1.0 + radial_distortion / mag};
    } else if (original_radius > r_max_ + eps_) {
      return std::nullopt;
    }
  }

  // We are within r_min and r_max and not at a boundary
  return std::optional{original_radius / mag};
}

std::array<double, 3> SphereTransition::gradient(
    const std::array<double, 3>& source_coords) const {
  const double mag = magnitude(source_coords);

  // Short circuit for the interior
  if (interior_) {
    return std::array{0.0, 0.0, 0.0};
  }

  if (UNLIKELY(equal_within_roundoff(mag, 0.0))) {
    ERROR("Trying to divide by a point "
          << source_coords
          << " with radius zero in SphereTransition gradient.");
  }

  if (UNLIKELY(mag < r_min_ - eps_)) {
    ERROR("SphereTransition gradient coord "
          << source_coords << " with radius " << mag << " is within r_min of "
          << r_min_
          << ", but the class was not constructed for the interior of the "
             "sphere.");
  } else if (UNLIKELY(mag > r_max_ + eps_)) {
    ERROR("SphereTransition gradient coord " << source_coords << " with radius "
                                             << mag << "is beyond r_max of "
                                             << r_max_ << ".");
  }

  // We can call the operator() and be sure it won't error because we did the
  // checks here as well.
  return source_coords * (a_ / square(mag) - (*this)(source_coords, {2}));
}
std::array<DataVector, 3> SphereTransition::gradient(
    const std::array<DataVector, 3>& source_coords) const {
  auto result = source_coords;
  for (size_t i = 0; i < source_coords[0].size(); i++) {
    auto double_result = gradient(std::array{
        source_coords[0][i], source_coords[1][i], source_coords[2][i]});
    for (size_t j = 0; j < 3; j++) {
      gsl::at(result, j)[i] = gsl::at(double_result, j);
    }
  }
  return result;
}

bool SphereTransition::operator==(
    const ShapeMapTransitionFunction& other) const {
  if (dynamic_cast<const SphereTransition*>(&other) == nullptr) {
    return false;
  }
  const auto& derived = dynamic_cast<const SphereTransition&>(other);
  // no need to check `a_` or `b_` as they are uniquely determined by `r_min_`
  // and `r_max_`.
  return this->r_min_ == derived.r_min_ and this->r_max_ == derived.r_max_ and
         this->interior_ == derived.interior_;
}

bool SphereTransition::operator!=(
    const ShapeMapTransitionFunction& other) const {
  return not(*this == other);
}

void SphereTransition::pup(PUP::er& p) {
  ShapeMapTransitionFunction::pup(p);
  size_t version = 1;
  p | version;
  // Remember to increment the version number when making changes to this
  // function. Retain support for unpacking data written by previous versions
  // whenever possible. See `Domain` docs for details.
  if (version >= 0) {
    p | r_min_;
    p | r_max_;
    p | a_;
    p | b_;
    if (version >= 1) {
      p | interior_;
    } else {
      interior_ = false;
    }
  } else if (p.isUnpacking()) {
    interior_ = false;
  }
}

SphereTransition::SphereTransition(CkMigrateMessage* const msg)
    : ShapeMapTransitionFunction(msg) {}

PUP::able::PUP_ID SphereTransition::my_PUP_ID = 0;  // NOLINT

}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
