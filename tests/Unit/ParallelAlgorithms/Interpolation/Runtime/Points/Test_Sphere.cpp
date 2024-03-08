// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/Runtime/PointsHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/Sphere.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/Spherepack.hpp"

namespace {
template <typename Generator>
void test_sphere(const gsl::not_null<Generator*> gen,
                 const size_t number_of_spheres,
                 const intrp2::AngularOrdering angular_ordering) {
  std::uniform_real_distribution<double> dist{1.2, 4.5};
  std::set<double> radii_set{};
  for (size_t i = 0; i < number_of_spheres; i++) {
    double radius = dist(*gen);
    while (radii_set.contains(radius)) {
      radius = dist(*gen);
    }
    radii_set.insert(radius);
  }
  const size_t l_max = 18;
  const std::array<double, 3> center = {{0.05, 0.06, 0.07}};

  CAPTURE(l_max);
  CAPTURE(center);
  CAPTURE(radii_set);
  CAPTURE(number_of_spheres);
  CAPTURE(angular_ordering);

  // Options for Sphere
  std::string radii_str;
  std::stringstream ss;
  ss << std::setprecision(std::numeric_limits<double>::max_digits10);
  if (number_of_spheres == 1) {
    // Test the double variant
    ss << *radii_set.begin();
  } else {
    // Test the vector variant
    auto it = radii_set.begin();
    ss << "[" << *it;
    it++;
    for (; it != radii_set.end(); it++) {
      ss << "," << *it;
    }
    ss << "]";
  }
  radii_str = ss.str();

  const std::string option_string =
      "Center: [0.05, 0.06, 0.07]\n"
      "Radius: " +
      radii_str +
      "\n"
      "LMax: 18\n"
      "AngularOrdering: " +
      std::string(MakeString{} << angular_ordering);

  // How many points are supposed to be in a Strahlkorper,
  // reproduced here by hand for the test.
  const size_t n_theta = l_max + 1;
  const size_t n_phi = 2 * l_max + 1;

  tnsr::I<DataVector, 3, Frame::NoFrame> expected_points{number_of_spheres *
                                                         n_theta * n_phi};

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
    if (angular_ordering == intrp2::AngularOrdering::Strahlkorper) {
      for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
        const double phi = two_pi_over_n_phi * i_phi;
        for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
          const double theta = theta_points[i_theta];
          expected_points.get(0)[s] =
              radius * sin(theta) * cos(phi) + center[0];
          expected_points.get(1)[s] =
              radius * sin(theta) * sin(phi) + center[1],
          expected_points.get(2)[s] = radius * cos(theta) + center[2];
          ++s;
        }
      }
    } else {
      for (size_t i_theta = 0; i_theta < n_theta; ++i_theta) {
        for (size_t i_phi = 0; i_phi < n_phi; ++i_phi) {
          const double phi = two_pi_over_n_phi * i_phi;
          const double theta = theta_points[i_theta];
          expected_points.get(0)[s] =
              radius * sin(theta) * cos(phi) + center[0];
          expected_points.get(1)[s] =
              radius * sin(theta) * sin(phi) + center[1],
          expected_points.get(2)[s] = radius * cos(theta) + center[2];
          ++s;
        }
      }
    }
  }

  const auto sphere =
      RuntimeIntrpTestHelpers::test_points<intrp2::points::Sphere>(
          option_string, expected_points);

  CHECK(sphere.l_max() == l_max);
  CHECK(sphere.center() == center);
  CHECK(sphere.radii() == radii_set);
  CHECK(sphere.angular_ordering() == angular_ordering);
}

// TODO: Still need to figure out the compute tag stuff for the Sphere points
// void test_operator() {
//   // We do this so we know exactly where some points are
//   const size_t l_max = 2;
//   const std::array<double, 3> center{0.0, 0.0, 0.0};
//   const std::set<double> radii_set{1.0, 2.0};
//   const intrp2::AngularOrdering angular_ordering =
//       intrp2::AngularOrdering::Strahlkorper;
//   const intrp2::points::Sphere sphere{
//       l_max, center, std::vector<double>{radii_set.begin(), radii_set.end()},
//       angular_ordering};
// }

void test_sphere_errors() {
  CHECK_THROWS_WITH(
      ([]() {
        const auto created_opts =
            TestHelpers::test_creation<intrp2::points::Sphere>(
                "Center: [0.05, 0.06, 0.07]\n"
                "Radius: [1.0, 1.0]\n"
                "LMax: 18\n"
                "AngularOrdering: Cce");
      })(),
      Catch::Matchers::ContainsSubstring(
          "into radii for Sphere interpolation target. It already "
          "exists. Existing radii are"));
  CHECK_THROWS_WITH(
      ([]() {
        const auto created_opts =
            TestHelpers::test_creation<intrp2::points::Sphere>(
                "Center: [0.05, 0.06, 0.07]\n"
                "Radius: [-1.0]\n"
                "LMax: 18\n"
                "AngularOrdering: Cce");
      })(),
      Catch::Matchers::ContainsSubstring("Radius must be positive"));
  CHECK_THROWS_WITH(
      ([]() {
        const auto created_opts =
            TestHelpers::test_creation<intrp2::points::Sphere>(
                "Center: [0.05, 0.06, 0.07]\n"
                "Radius: -1.0\n"
                "LMax: 18\n"
                "AngularOrdering: Cce");
      })(),
      Catch::Matchers::ContainsSubstring("Radius must be positive"));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ParallelAlgorithms.Interpolation.Runtime.Points.Sphere",
                  "[Unit]") {
  MAKE_GENERATOR(gen);
  test_sphere_errors();
  for (const auto& [num_spheres, angular_ordering] :
       cartesian_product(std::array{1_st, 2_st, 3_st},
                         std::array{intrp2::AngularOrdering::Cce,
                                    intrp2::AngularOrdering::Strahlkorper})) {
    test_sphere(make_not_null(&gen), num_spheres, angular_ordering);
  }
}
