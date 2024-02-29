// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Runtime/Points/Sphere.hpp"

#include <array>
#include <cstddef>
#include <set>
#include <string>
#include <variant>

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Transpose.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/StrahlkorperFunctions.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"
#include "Utilities/StdHelpers.hpp"

namespace intrp2::points {
namespace {
struct SphereVisitor {
  const Options::Context& context;

  std::set<double> operator()(const double radius) {
    positive_radius(radius);
    return std::set<double>{radius};
  }

  std::set<double> operator()(const std::vector<double>& radii) {
    std::set<double> result;
    for (const double radius : radii) {
      if (result.count(radius) != 0) {
        using ::operator<<;
        PARSE_ERROR(context,
                    "Cannot insert radius "
                        << radius
                        << " into radii for Sphere interpolation target. It "
                           "already exists. Existing radii are "
                        << result);
      }
      positive_radius(radius);
      result.emplace(radius);
    }
    return result;
  }

 private:
  void positive_radius(const double radius) {
    if (radius <= 0) {
      PARSE_ERROR(context, "Radius must be positive, not " << radius);
    }
  }
};
}  // namespace

Sphere::Sphere(const size_t l_max, const std::array<double, 3>& center,
               const typename Radius::type& input_radii,
               const intrp2::AngularOrdering angular_ordering,
               const Options::Context& context)
    : l_max_(l_max),
      center_(center),
      radii_(std::visit(SphereVisitor{context}, input_radii)),
      angular_ordering_(angular_ordering) {
  // Total number of points is number of points for one sphere times the
  // number of spheres we use.
  const size_t num_points_single_sphere = (l_max_ + 1) * (2 * l_max_ + 1);
  const size_t num_points = radii_.size() * num_points_single_sphere;

  set_number_of_grid_points(make_not_null(&points_), num_points);

  size_t index = 0;
  for (const double radius : radii_) {
    ylm::Strahlkorper<Frame::NoFrame> strahlkorper(
        l_max_, l_max_, DataVector{(l_max_ + 1) * (2 * l_max_ + 1), radius},
        center_);

    tnsr::I<DataVector, 3, Frame::NoFrame> single_sphere_coords =
        ylm::cartesian_coords(strahlkorper, ylm::radius(strahlkorper),
                              ylm::rhat(ylm::theta_phi(strahlkorper)));

    // If the angular ordering is Strahlkorper then we don't have to do
    // anything to the coords because they are already in the right order
    if (angular_ordering_ == intrp2::AngularOrdering::Cce) {
      const auto physical_extents =
          strahlkorper.ylm_spherepack().physical_extents();

      for (size_t i = 0; i < 3; ++i) {
        // Can't use not_null version because points get overwritten before they
        // are used
        single_sphere_coords.get(i) =
            transpose(single_sphere_coords.get(i), physical_extents[0],
                      physical_extents[1]);
      }
    }

    const size_t tmp_index = index;
    for (size_t i = 0; i < 3; ++i) {
      for (size_t local_index = 0;
           local_index < single_sphere_coords.get(i).size();
           local_index++, index++) {
        points_.get(i)[index] = single_sphere_coords.get(i)[local_index];
      }
      if (i != 2) {
        index = tmp_index;
      }
    }
  }

  // If this fails, there is a bug. Can't really test it
  // LCOV_EXCL_START
  ASSERT(index == points_.get(0).size(),
         "Didn't initialize points of Sphere target correctly. index = "
             << index << " all_coords.size() = " << points_.get(0).size());
  // LCOV_EXCL_STOP
}

size_t Sphere::l_max() const { return l_max_; }
const std::array<double, 3>& Sphere::center() const { return center_; }
const std::set<double>& Sphere::radii() const { return radii_; }
intrp2::AngularOrdering Sphere::angular_ordering() const {
  return angular_ordering_;
}

void Sphere::pup(PUP::er& p) {
  p | l_max_;
  p | center_;
  p | radii_;
  p | angular_ordering_;
  p | points_;
}

const tnsr::I<DataVector, 3, Frame::NoFrame>& Sphere::target_points_no_frame()
    const {
  return points_;
}

bool operator==(const Sphere& lhs, const Sphere& rhs) {
  return lhs.l_max_ == rhs.l_max_ and lhs.center_ == rhs.center_ and
         lhs.radii_ == rhs.radii_ and
         lhs.angular_ordering_ == rhs.angular_ordering_ and
         lhs.points_ == rhs.points_;
}

bool operator!=(const Sphere& lhs, const Sphere& rhs) {
  return not(lhs == rhs);
}
}  // namespace intrp2::points
