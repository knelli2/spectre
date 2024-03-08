// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/Runtime/PointsHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/KerrHorizon.hpp"
#include "Utilities/Spherepack.hpp"

namespace {
void test_kerr_horizon(const intrp2::AngularOrdering angular_ordering) {
  std::uniform_real_distribution<double> dist{1.2, 4.5};
  const size_t l_max = 18;
  const double mass = 1.8;
  const std::array<double, 3> center{0.05, 0.06, 0.07};
  const std::array<double, 3> dimensionless_spin{0.2, 0.3, 0.4};

  CAPTURE(l_max);
  CAPTURE(center);
  CAPTURE(dimensionless_spin);
  CAPTURE(angular_ordering);

  const std::string option_string =
      "Center: [0.05, 0.06, 0.07]\n"
      "DimensionlessSpin: [0.2, 0.3, 0.4]\n"
      "Mass: 1.8\n"
      "LMax: 18\n"
      "AngularOrdering: " +
      std::string(MakeString{} << angular_ordering);

  // How many points are supposed to be in a Strahlkorper,
  // reproduced here by hand for the test.
  const size_t n_theta = l_max + 1;
  const size_t n_phi = 2 * l_max + 1;

  tnsr::I<DataVector, 3, Frame::NoFrame> expected_points{n_theta * n_phi};

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
  if (angular_ordering == intrp2::AngularOrdering::Strahlkorper) {
    for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
      const double phi = two_pi_over_n_phi * i_phi;
      for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
        const double theta = theta_points[i_theta];
        const double r = radius(theta, phi);
        expected_points.get(0)[s] = r * sin(theta) * cos(phi) + center[0];
        expected_points.get(1)[s] = r * sin(theta) * sin(phi) + center[1],
        expected_points.get(2)[s] = r * cos(theta) + center[2];
        ++s;
      }
    }
  } else {
    for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
      for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
        const double phi = two_pi_over_n_phi * i_phi;
        const double theta = theta_points[i_theta];
        const double r = radius(theta, phi);
        expected_points.get(0)[s] = r * sin(theta) * cos(phi) + center[0];
        expected_points.get(1)[s] = r * sin(theta) * sin(phi) + center[1],
        expected_points.get(2)[s] = r * cos(theta) + center[2];
        ++s;
      }
    }
  }

  const auto kerr_horizon =
      RuntimeIntrpTestHelpers::test_points<intrp2::points::KerrHorizon>(
          option_string, expected_points);

  CHECK(kerr_horizon.l_max() == l_max);
  CHECK(kerr_horizon.center() == center);
  CHECK(kerr_horizon.mass() == mass);
  CHECK(kerr_horizon.dimensionless_spin() == dimensionless_spin);
  CHECK(kerr_horizon.angular_ordering() == angular_ordering);
}

void test_kerr_horizon_errors() {
  CHECK_THROWS_WITH(
      ([]() {
        const auto created_opts =
            TestHelpers::test_creation<intrp2::points::KerrHorizon>(
                "Center: [0.05, 0.06, 0.07]\n"
                "DimensionlessSpin: [0.2, 0.3, 0.4]\n"
                "Mass: -50.0\n"
                "LMax: 18\n"
                "AngularOrdering: Cce");
      })(),
      Catch::Matchers::ContainsSubstring("KerrHorizon expects mass > 0"));
  CHECK_THROWS_WITH(
      ([]() {
        const auto created_opts =
            TestHelpers::test_creation<intrp2::points::KerrHorizon>(
                "Center: [0.05, 0.06, 0.07]\n"
                "DimensionlessSpin: [1.2, 0.3, 0.4]\n"
                "Mass: 2.0\n"
                "LMax: 18\n"
                "AngularOrdering: Cce");
      })(),
      Catch::Matchers::ContainsSubstring(
          "KerrHorizon expects |dimensionless_spin|<=1"));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.Points.KerrHorizon",
    "[Unit]") {
  test_kerr_horizon_errors();
  for (const auto& angular_ordering :
       std::array{intrp2::AngularOrdering::Cce,
                  intrp2::AngularOrdering::Strahlkorper}) {
    test_kerr_horizon(angular_ordering);
  }
}
