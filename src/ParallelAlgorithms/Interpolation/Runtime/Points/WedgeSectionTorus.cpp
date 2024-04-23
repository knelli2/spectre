// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Runtime/Points/WedgeSectionTorus.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace intrp2::points {
WedgeSectionTorus::WedgeSectionTorus(
    const std::array<double, 2>& radial_bounds,
    const std::array<double, 2>& theta_bounds,
    const std::array<size_t, 3>& number_of_grid_points,
    const bool use_uniform_radial_grid, const bool use_uniform_theta_grid,
    const Options::Context& context)
    : radial_bounds_(radial_bounds),
      theta_bounds_(theta_bounds),
      number_of_grid_points_(number_of_grid_points),
      use_uniform_radial_grid_(use_uniform_radial_grid),
      use_uniform_theta_grid_(use_uniform_theta_grid) {
  const double min_radius = gsl::at(radial_bounds_, 0);
  const double max_radius = gsl::at(radial_bounds_, 1);
  const double min_theta = gsl::at(theta_bounds_, 0);
  const double max_theta = gsl::at(theta_bounds_, 1);
  const size_t number_of_phi_points = gsl::at(number_of_grid_points_, 0);
  const size_t number_of_theta_points = gsl::at(number_of_grid_points_, 1);
  const size_t number_of_radial_points = gsl::at(number_of_grid_points_, 2);

  if (min_radius < 0.0 or max_radius < 0.0) {
    PARSE_ERROR(context, "Both radial bounds must be non-negative.");
  }
  if (min_theta < 0.0 or max_theta < 0.0) {
    PARSE_ERROR(context, "Both theta bounds must be non-negative.");
  }
  if (min_theta > M_PI or max_theta > M_PI) {
    PARSE_ERROR(context, "Both theta bounds cannot be larger than pi.");
  }
  if (number_of_phi_points < 1) {
    PARSE_ERROR(context, "Number of phi points must be larger than 0.");
  }
  if (number_of_theta_points < 2) {
    PARSE_ERROR(context, "Number of theta points must be larger than 1.");
  }
  if (number_of_radial_points < 2) {
    PARSE_ERROR(context, "Number of radial points must be larger than 1.");
  }
  if (min_radius >= max_radius) {
    PARSE_ERROR(context, "WedgeSectionTorus expects min_radius < max_radius.");
  }
  if (min_theta >= max_theta) {
    PARSE_ERROR(context, "WedgeSectionTorus expects min_theta < max_theta.");
  }

  // Compute locations of constant r/theta/phi surfaces
  const DataVector radii_1d = [this, &number_of_radial_points, &max_radius,
                               &min_radius]() {
    DataVector result(number_of_radial_points);
    if (use_uniform_radial_grid_) {
      // uniform point distribution
      const double coefficient =
          (max_radius - min_radius) / (number_of_radial_points - 1.0);
      for (size_t r = 0; r < number_of_radial_points; ++r) {
        result[r] = min_radius + coefficient * r;
      }
    } else {
      // radial Legendre-Gauss-Lobatto distribution via linear rescaling
      const double mean = 0.5 * (max_radius + min_radius);
      const double diff = 0.5 * (max_radius - min_radius);
      result =
          mean +
          diff *
              Spectral::collocation_points<Spectral::Basis::Legendre,
                                           Spectral::Quadrature::GaussLobatto>(
                  number_of_radial_points);
    }
    return result;
  }();
  const DataVector thetas_1d = [this, &number_of_theta_points, &max_theta,
                                &min_theta]() {
    DataVector result(number_of_theta_points);
    if (use_uniform_theta_grid_) {
      // uniform point distribution
      const double coefficient =
          (max_theta - min_theta) / (number_of_theta_points - 1.0);
      for (size_t theta = 0; theta < number_of_theta_points; ++theta) {
        result[theta] = min_theta + coefficient * theta;
      }
    } else {
      // Legendre-Gauss-Lobatto point distribution (in theta)
      const double mean = 0.5 * (max_theta + min_theta);
      const double diff = 0.5 * (max_theta - min_theta);
      result =
          mean +
          diff *
              Spectral::collocation_points<Spectral::Basis::Legendre,
                                           Spectral::Quadrature::GaussLobatto>(
                  number_of_theta_points);
    }
    return result;
  }();
  const DataVector phis_1d = [&number_of_phi_points]() {
    DataVector result(number_of_phi_points);
    // We do NOT want a grid point at phi = 2pi, as this would duplicate the
    // phi = 0 data. So, divide by number_of_phi_points rather than (n-1) as
    // elsewhere.
    const double coefficient = 2.0 * M_PI / number_of_phi_points;
    for (size_t phi = 0; phi < number_of_phi_points; ++phi) {
      result[phi] = coefficient * phi;
    }
    return result;
  }();

  // Take tensor product to get full 3D r/theta/phi points
  const size_t num_total =
      number_of_radial_points * number_of_theta_points * number_of_phi_points;
  DataVector radii(num_total);
  DataVector thetas(num_total);
  DataVector phis(num_total);
  for (size_t phi = 0; phi < number_of_phi_points; ++phi) {
    for (size_t theta = 0; theta < number_of_theta_points; ++theta) {
      for (size_t r = 0; r < number_of_radial_points; ++r) {
        const size_t i = r + theta * number_of_radial_points +
                         phi * number_of_theta_points * number_of_radial_points;
        radii[i] = radii_1d[r];
        thetas[i] = thetas_1d[theta];
        phis[i] = phis_1d[phi];
      }
    }
  }

  // Compute x/y/z coordinates
  // Note: theta measured from +z axis, phi measured from +x axis
  set_number_of_grid_points(make_not_null(&points_), num_total);
  get<0>(points_) = radii * sin(thetas) * cos(phis);
  get<1>(points_) = radii * sin(thetas) * sin(phis);
  get<2>(points_) = radii * cos(thetas);
}

const std::array<double, 2>& WedgeSectionTorus::radial_bounds() const {
  return radial_bounds_;
}
const std::array<double, 2>& WedgeSectionTorus::theta_bounds() const {
  return theta_bounds_;
}
const std::array<size_t, 3>& WedgeSectionTorus::number_of_grid_points() const {
  return number_of_grid_points_;
}
bool WedgeSectionTorus::use_uniform_radial_grid() const {
  return use_uniform_radial_grid_;
}
bool WedgeSectionTorus::use_uniform_theta_grid() const {
  return use_uniform_theta_grid_;
}

void WedgeSectionTorus::pup(PUP::er& p) {
  p | radial_bounds_;
  p | theta_bounds_;
  p | number_of_grid_points_;
  p | use_uniform_radial_grid_;
  p | use_uniform_theta_grid_;
  p | points_;
}

const tnsr::I<DataVector, 3, Frame::NoFrame>&
WedgeSectionTorus::target_points_no_frame() const {
  return points_;
}

size_t WedgeSectionTorus::number_of_sets_of_points() const { return 1; }

bool operator==(const WedgeSectionTorus& lhs, const WedgeSectionTorus& rhs) {
  return lhs.radial_bounds_ == rhs.radial_bounds_ and
         lhs.theta_bounds_ == rhs.theta_bounds_ and
         lhs.number_of_grid_points_ == rhs.number_of_grid_points_ and
         lhs.use_uniform_radial_grid_ == rhs.use_uniform_radial_grid_ and
         lhs.use_uniform_theta_grid_ == rhs.use_uniform_theta_grid_ and
         lhs.points_ == rhs.points_;
}

bool operator!=(const WedgeSectionTorus& lhs, const WedgeSectionTorus& rhs) {
  return not(lhs == rhs);
}
}  // namespace intrp2::points
