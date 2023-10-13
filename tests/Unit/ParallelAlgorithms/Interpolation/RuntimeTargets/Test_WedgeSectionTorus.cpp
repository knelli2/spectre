// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <memory>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/RuntimeTargetHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Target.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/WedgeSectionTorus.hpp"

namespace intrp::Targets {
namespace {
void test_wedge_section_torus() {
  const std::string expected_name = "WedgeSectionTorus";
  const std::string option_string =
      "WedgeSectionTorus:\n"
      " RadialBounds: [0.1, 0.5]\n"
      " ThetaBounds: [1.8, 2.4]\n"
      " NumberOfGridPoints: [8, 4, 4]\n"
      " UniformRadialGrid: false\n"
      " UniformThetaGrid: false\n";

  const size_t num_radial = 4;
  const size_t num_theta = 4;
  const size_t num_phi = 8;
  const size_t num_total = num_radial * num_theta * num_phi;
  tnsr::I<DataVector, 3, Frame::NoFrame> expected_points_no_frame{num_total};
  for (size_t r = 0; r < num_radial; ++r) {
    const double radius =
        0.3 +
        0.2 * Spectral::collocation_points<Spectral::Basis::Legendre,
                                           Spectral::Quadrature::GaussLobatto>(
                  num_radial)[r];
    for (size_t t = 0; t < num_theta; ++t) {
      const double theta =
          2.1 +
          0.3 *
              Spectral::collocation_points<Spectral::Basis::Legendre,
                                           Spectral::Quadrature::GaussLobatto>(
                  num_theta)[t];
      for (size_t p = 0; p < num_phi; ++p) {
        const double phi = 2.0 * M_PI * p / num_phi;
        const size_t i = r + t * num_radial + p * num_theta * num_radial;
        get<0>(expected_points_no_frame)[i] = radius * sin(theta) * cos(phi);
        get<1>(expected_points_no_frame)[i] = radius * sin(theta) * sin(phi);
        get<2>(expected_points_no_frame)[i] = radius * cos(theta);
      }
    }
  }

  const auto target = TestHelpers::test_target_all_frames<WedgeSectionTorus>(
      option_string, expected_points_no_frame, expected_name);
  const WedgeSectionTorus& torus =
      dynamic_cast<const WedgeSectionTorus&>(*target.get());

  CHECK(torus.radial_bounds() == std::array{0.1, 0.5});
  CHECK(torus.theta_bounds() == std::array{1.8, 2.4});
  CHECK(torus.number_of_grid_points() == std::array{8_st, 4_st, 4_st});
  CHECK_FALSE(torus.use_uniform_radial_grid());
  CHECK_FALSE(torus.use_uniform_theta_grid());
}

void test_errors() {
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{-1.0, 2.0}, std::array{1.0, 2.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Both radial bounds must be non-negative."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, -2.0}, std::array{1.0, 2.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Both radial bounds must be non-negative."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{-1.0, 2.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds must be non-negative."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{1.0, -2.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds must be non-negative."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{3.15, 2.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds cannot be larger than pi."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{1.0, 3.15},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds cannot be larger than pi."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{1.0, 2.0},
                         std::array{0_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Number of phi points must be larger than 0."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{1.0, 2.0},
                         std::array{1_st, 1_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Number of theta points must be larger than 1."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{1.0, 2.0},
                         std::array{1_st, 2_st, 1_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "Number of radial points must be larger than 1."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{2.0, 1.0}, std::array{1.0, 2.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "WedgeSectionTorus expects min_radius < max_radius."));
  CHECK_THROWS_WITH(
      (WedgeSectionTorus{std::array{1.0, 2.0}, std::array{2.0, 1.0},
                         std::array{1_st, 2_st, 2_st}, false, false}),
      Catch::Matchers::ContainsSubstring(
          "WedgeSectionTorus expects min_theta < max_theta."));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.RuntimeTargets.WedgeSectionTorus",
    "[Unit]") {
  test_errors();
  test_wedge_section_torus();
}
}  // namespace intrp::Targets
