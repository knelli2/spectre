// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/Runtime/PointsHelpers.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/WedgeSectionTorus.hpp"

namespace {
void test_wedge_section_torus() {
  const std::string option_string{
      "RadialBounds: [1.2, 4.0]\n"
      "ThetaBounds: [1.1, 1.73]\n"
      "NumberOfGridPoints: [8, 4, 8]\n"
      "UniformRadialGrid: false\n"
      "UniformThetaGrid: false\n"};

  const size_t num_radial = 8;
  const size_t num_theta = 4;
  const size_t num_phi = 8;
  const double min_radius = 1.2;
  const double max_radius = 4.0;
  const double avg_radius = 0.5 * (max_radius + min_radius);
  const double diff_radius = 0.5 * (max_radius - min_radius);
  const double min_theta = 1.1;
  const double max_theta = 1.73;
  const double avg_theta = 0.5 * (max_theta + min_theta);
  const double diff_theta = 0.5 * (max_theta - min_theta);
  const size_t num_total = num_radial * num_theta * num_phi;
  tnsr::I<DataVector, 3, Frame::NoFrame> expected_points{num_total};
  for (size_t r = 0; r < num_radial; ++r) {
    const double radius =
        avg_radius +
        diff_radius *
            Spectral::collocation_points<Spectral::Basis::Legendre,
                                         Spectral::Quadrature::GaussLobatto>(
                num_radial)[r];
    for (size_t t = 0; t < num_theta; ++t) {
      const double theta =
          avg_theta +
          diff_theta *
              Spectral::collocation_points<Spectral::Basis::Legendre,
                                           Spectral::Quadrature::GaussLobatto>(
                  num_theta)[t];
      for (size_t p = 0; p < num_phi; ++p) {
        const double phi = 2.0 * M_PI * p / num_phi;
        const size_t i = r + t * num_radial + p * num_theta * num_radial;
        get<0>(expected_points)[i] = radius * sin(theta) * cos(phi);
        get<1>(expected_points)[i] = radius * sin(theta) * sin(phi);
        get<2>(expected_points)[i] = radius * cos(theta);
      }
    }
  }

  const auto torus =
      RuntimeIntrpTestHelpers::test_points<intrp2::points::WedgeSectionTorus>(
          option_string, expected_points);

  CHECK(torus.number_of_sets_of_points() == 1);
}

void test_wedge_section_torus_errors() {
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {-1.2, 4.0}, {1.1, 1.73}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Both radial bounds must be non-negative."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, -4.0}, {1.1, 1.73}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Both radial bounds must be non-negative."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {-1.1, 1.73}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds must be non-negative."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {1.1, -1.73}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds must be non-negative."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {4.0, 1.73}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds cannot be larger than pi."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {1.1, 4.0}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Both theta bounds cannot be larger than pi."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {1.1, 1.73}, {0_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Number of phi points must be larger than 0."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {1.1, 1.73}, {8_st, 0_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Number of theta points must be larger than 1."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {1.1, 1.73}, {8_st, 4_st, 0_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "Number of radial points must be larger than 1."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {4.0, 1.2}, {1.1, 1.73}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "WedgeSectionTorus expects min_radius < max_radius."));
  CHECK_THROWS_WITH(
      ([]() {
        const intrp2::points::WedgeSectionTorus torus{
            {1.2, 4.0}, {1.73, 1.1}, {8_st, 4_st, 8_st}, false, true};
        (void)torus;
      })(),
      Catch::Matchers::ContainsSubstring(
          "WedgeSectionTorus expects min_theta < max_theta."));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.Points.WedgeSectionTorus",
    "[Unit]") {
  test_wedge_section_torus_errors();
  test_wedge_section_torus();
}
