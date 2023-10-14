// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/RuntimeTargetHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/KerrHorizon.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Spherepack.hpp"

namespace intrp::Targets {
namespace {
void test_kerr_horizon(const intrp::AngularOrdering angular_ordering) {
  const std::string expected_name = "KerrHorizon";
  const std::string option_string =
      "KerrHorizon:\n"
      " LMax: 4\n"
      " Center: [0.1, 0.3, 0.4]\n"
      " Mass: 0.3\n"
      " DimensionlessSpin: [0.0, 0.0, 0.2]\n"
      " AngularOrdering: " +
      get_output(angular_ordering) + "\n";

  // How many points are supposed to be in a Strahlkorper,
  // reproduced here by hand for the test.
  const size_t l_max = 4;
  const std::array<double, 3> center{0.1, 0.3, 0.4};
  const double mass = 0.3;
  const std::array<double, 3> dimensionless_spin{0.0, 0.0, 0.2};
  const size_t n_theta = l_max + 1;
  const size_t n_phi = 2 * l_max + 1;

  tnsr::I<DataVector, 3, Frame::NoFrame> expected_points_no_frame{n_theta *
                                                                  n_phi};

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

  // Radius as function of theta, phi
  const auto radius = [&mass, &dimensionless_spin](const double theta,
                                                   const double phi) {
    // Recoding kerr_horizon_radius in a different way for the test.
    const std::array<double, 3> spin_a = {{mass * dimensionless_spin[0],
                                           mass * dimensionless_spin[1],
                                           mass * dimensionless_spin[2]}};
    const double spin_a_squared =
        square(spin_a[0]) + square(spin_a[1]) + square(spin_a[2]);
    const double a_dot_xhat_squared =
        square(spin_a[0] * sin(theta) * cos(phi) +
               spin_a[1] * sin(theta) * sin(phi) + spin_a[2] * cos(theta));
    const double r_boyer_lindquist_squared =
        square(mass + sqrt(square(mass) - spin_a_squared));
    return sqrt((r_boyer_lindquist_squared + spin_a_squared) /
                (1.0 + a_dot_xhat_squared / r_boyer_lindquist_squared));
  };

  const double two_pi_over_n_phi = 2.0 * M_PI / n_phi;
  size_t s = 0;
  if (angular_ordering == intrp::AngularOrdering::Strahlkorper) {
    for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
      const double phi = two_pi_over_n_phi * i_phi;
      for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
        const double theta = theta_points[i_theta];
        const double r = radius(theta, phi);
        expected_points_no_frame.get(0)[s] =
            r * sin(theta) * cos(phi) + center[0];
        expected_points_no_frame.get(1)[s] =
            r * sin(theta) * sin(phi) + center[1],
        expected_points_no_frame.get(2)[s] = r * cos(theta) + center[2];
        ++s;
      }
    }
  } else {
    for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
      for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
        const double phi = two_pi_over_n_phi * i_phi;
        const double theta = theta_points[i_theta];
        const double r = radius(theta, phi);
        expected_points_no_frame.get(0)[s] =
            r * sin(theta) * cos(phi) + center[0];
        expected_points_no_frame.get(1)[s] =
            r * sin(theta) * sin(phi) + center[1],
        expected_points_no_frame.get(2)[s] = r * cos(theta) + center[2];
        ++s;
      }
    }
  }

  const auto target = TestHelpers::test_target_all_frames<KerrHorizon>(
      option_string, expected_points_no_frame, expected_name);
  const KerrHorizon& kerr_horizon =
      dynamic_cast<const KerrHorizon&>(*target.get());

  CHECK(kerr_horizon.l_max() == l_max);
  CHECK(kerr_horizon.center() == center);
  CHECK(kerr_horizon.mass() == mass);
  CHECK(kerr_horizon.dimensionless_spin() == dimensionless_spin);
  CHECK(kerr_horizon.angular_ordering() == angular_ordering);
}

void test_errors() {
  CHECK_THROWS_WITH(
      (KerrHorizon{2_st, std::array{0.0, 0.0, 0.0}, -1.0,
                   std::array{0.0, 0.0, 0.0}, AngularOrdering::Strahlkorper}),
      Catch::Matchers::ContainsSubstring(
          "KerrHorizon expects mass > 0, not -1"));
  CHECK_THROWS_WITH(
      (KerrHorizon{2_st, std::array{0.0, 0.0, 0.0}, 1.0,
                   std::array{0.0, 0.0, 2.0}, AngularOrdering::Strahlkorper}),
      Catch::Matchers::ContainsSubstring(
          "KerrHorizon expects |dimensionless_spin|<=1, not 2"));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.RuntimeTargets.KerrHorizon",
    "[Unit]") {
  test_errors();
  test_kerr_horizon(intrp::AngularOrdering::Strahlkorper);
  test_kerr_horizon(intrp::AngularOrdering::Cce);
}
}  // namespace intrp::Targets
