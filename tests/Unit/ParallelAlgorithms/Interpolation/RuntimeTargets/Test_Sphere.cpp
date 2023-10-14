// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/RuntimeTargetHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Sphere.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Spherepack.hpp"

namespace intrp::Targets {
namespace {
void test_sphere(const intrp::AngularOrdering angular_ordering) {
  const std::string expected_name = "Sphere";
  const std::string option_string =
      "Sphere:\n"
      " LMax: 4\n"
      " Center: [-0.1, -0.2, -0.3]\n"
      " Radius: [1.2, 2.3]\n"
      " AngularOrdering: " +
      get_output(angular_ordering) + "\n";

  const size_t number_of_spheres = 2;
  const size_t l_max = 4;
  const std::array<double, 3> center{-0.1, -0.2, -0.3};
  // How many points are supposed to be in a Strahlkorper,
  // reproduced here by hand for the test.
  const size_t n_theta = l_max + 1;
  const size_t n_phi = 2 * l_max + 1;

  // Have to turn this into a set to guarantee ordering
  const std::set<double> radii_set{1.2, 2.3};

  tnsr::I<DataVector, 3, Frame::NoFrame> expected_points_no_frame(
      number_of_spheres * n_theta * n_phi);

  size_t s = 0;
  for (const double radius : radii_set) {
    // The theta points of a Strahlkorper are Gauss-Legendre points.
    const std::vector<double> theta_points = []() {
      std::vector<double> thetas(n_theta);
      std::vector<double> work(n_theta + 1);
      std::vector<double> unused_weights(n_theta);
      int err = 0;
      gaqd_(static_cast<int>(n_theta), thetas.data(), unused_weights.data(),
            work.data(), static_cast<int>(n_theta + 1), &err);
      return thetas;
    }();

    const double two_pi_over_n_phi = 2.0 * M_PI / n_phi;
    if (angular_ordering == intrp::AngularOrdering::Strahlkorper) {
      for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
        const double phi = two_pi_over_n_phi * i_phi;
        for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
          const double theta = theta_points[i_theta];
          expected_points_no_frame.get(0)[s] =
              radius * sin(theta) * cos(phi) + center[0];
          expected_points_no_frame.get(1)[s] =
              radius * sin(theta) * sin(phi) + center[1],
          expected_points_no_frame.get(2)[s] = radius * cos(theta) + center[2];
          ++s;
        }
      }
    } else {
      for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
        for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
          const double phi = two_pi_over_n_phi * i_phi;
          const double theta = theta_points[i_theta];
          expected_points_no_frame.get(0)[s] =
              radius * sin(theta) * cos(phi) + center[0];
          expected_points_no_frame.get(1)[s] =
              radius * sin(theta) * sin(phi) + center[1],
          expected_points_no_frame.get(2)[s] = radius * cos(theta) + center[2];
          ++s;
        }
      }
    }
  }

  const auto target = TestHelpers::test_target_all_frames<Sphere>(
      option_string, expected_points_no_frame, expected_name);
  const Sphere& sphere = dynamic_cast<const Sphere&>(*target.get());

  CHECK(sphere.l_max() == l_max);
  CHECK(sphere.center() == center);
  CHECK(sphere.radii() == radii_set);
  CHECK(sphere.angular_ordering() == angular_ordering);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ParallelAlgorithms.Interpolation.RuntimeTargets.Sphere",
                  "[Unit]") {
  test_sphere(intrp::AngularOrdering::Strahlkorper);
  test_sphere(intrp::AngularOrdering::Cce);
}
}  // namespace intrp::Targets
