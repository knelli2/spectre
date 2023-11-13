// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "Domain/CoordinateMaps/TimeDependent/Shape.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Wedge.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/StdArrayHelpers.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {
namespace {
std::array<double, 3> sph_to_cart(const double radius, const double theta,
                                  const double phi) {
  return std::array{radius * sin(theta) * cos(phi),
                    radius * sin(theta) * sin(phi), radius * cos(theta)};
}

struct Expected {
  double inner_radius{};
  double outer_radius{};
  double inner_sphericity{};
  double outer_sphericity{};
  Falloff falloff{};
  std::vector<std::array<double, 3>> points_to_check{};
  std::vector<double> distortion{};
  std::vector<double> expected_transition{};
  std::vector<std::array<double, 3>> expected_mapped_points{};
  std::vector<std::array<double, 3>> expected_gradient{};

  void clear() {
    inner_sphericity = std::numeric_limits<double>::signaling_NaN();
    outer_sphericity = std::numeric_limits<double>::signaling_NaN();
    points_to_check.clear();
    distortion.clear();
    expected_transition.clear();
    expected_mapped_points.clear();
    expected_gradient.clear();
  }
};

template <typename T>
std::array<T, 2> cartesian_to_spherical(const std::array<T, 3>& cartesian) {
  const auto& [x, y, z] = cartesian;
  return {atan2(hypot(x, y), z), atan2(y, x)};
}

double compute_distortion(const gsl::not_null<ylm::Spherepack*> ylm,
                          const std::array<double, 3> coords,
                          const DataVector& coefs) {
  auto theta_phis = cartesian_to_spherical(coords);
  ylm->set_up_interpolation_info(theta_phis);

  return ylm->interpolate_from_coefs(coefs, theta_phis);
}

template <typename Generator>
void check_specific_points(const gsl::not_null<Generator*> generator) {
  INFO("Test specific points using shape map");
  std::uniform_real_distribution<double> coef_dist{-0.01, 0.01};

  const double initial_time = 1.0;
  const size_t l_max = 4;
  const size_t num_coefs = 2 * square(l_max + 1);
  ylm::Spherepack ylm{l_max, l_max};
  const std::array center{0.0, 0.0, 0.0};
  const std::string fot_name{"TheBean"};

  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      functions_of_time{};
  auto coefs = make_with_random_values<DataVector>(
      generator, make_not_null(&coef_dist), DataVector{num_coefs});
  functions_of_time[fot_name] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          initial_time, std::array{coefs},
          std::numeric_limits<double>::infinity());

  const auto test_points =
      [&](const Expected& local_expected, const double time,
          const std::unordered_map<
              std::string,
              std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
              functions_of_time) {
        std::unique_ptr<ShapeMapTransitionFunction> wedge =
            std::make_unique<Wedge>(
                local_expected.inner_radius, local_expected.outer_radius,
                local_expected.inner_sphericity,
                local_expected.outer_sphericity, local_expected.falloff);

        TimeDependent::Shape shape{center, l_max, l_max, wedge->get_clone(),
                                   fot_name};

        const size_t size = local_expected.points_to_check.size();
        REQUIRE(local_expected.distortion.size() == size);
        REQUIRE(local_expected.expected_transition.size() == size);
        REQUIRE(local_expected.expected_mapped_points.size() == size);
        REQUIRE(local_expected.expected_gradient.size() == size);

        for (size_t i = 0; i < local_expected.points_to_check.size(); i++) {
          const auto& point = local_expected.points_to_check[i];
          CAPTURE(i);
          CAPTURE(point);
          // Wedge functions
          CHECK(wedge->operator()(point) ==
                approx(local_expected.expected_transition[i]));
          const auto r_over_rtil_opt = wedge->original_radius_over_radius(
              local_expected.expected_mapped_points[i],
              local_expected.distortion[i]);
          CHECK(r_over_rtil_opt.has_value());
          const double r_over_rtil = r_over_rtil_opt.value();
          CHECK(magnitude(point) /
                    magnitude(local_expected.expected_mapped_points[i]) ==
                approx(r_over_rtil));
          const std::array gradient = wedge->gradient(point);
          CHECK_ITERABLE_APPROX(local_expected.expected_gradient[i], gradient);

          // Shape map. Not calculating jacobian
          const std::array mapped_point = shape(point, time, functions_of_time);
          CHECK_ITERABLE_APPROX(local_expected.expected_mapped_points[i],
                                mapped_point);
          const auto inverse_point =
              shape.inverse(mapped_point, time, functions_of_time);
          CHECK(inverse_point.has_value());
          CHECK_ITERABLE_APPROX(point, inverse_point.value());
        }
      };

  // We'll keep the radii the same since they don't really matter
  const double inner_radius = 0.5;
  const double outer_radius = 10.0;
  const double one_over_distance_difference =
      1.0 / (outer_radius - inner_radius);
  Expected expected{inner_radius, outer_radius};

  const auto compute_transition =
      [](const double inner_distance, const double outer_distance,
         const double radius, const Falloff falloff) {
        const double one_over_denom = 1.0 / (outer_distance - inner_distance);
        return falloff == Falloff::Linear
                   ? (outer_distance - radius) * one_over_denom
                   : -inner_distance * one_over_denom *
                         (1.0 - outer_distance / radius);
      };

  // We manually go through case by case since we know these points
  {
    INFO("Both spheres.");
    // Gradient is simple here
    const auto compute_gradient = [&](const std::array<double, 3>& point,
                                      const Falloff falloff) {
      const double radius = magnitude(point);
      const double extra_factor =
          falloff == Falloff::Linear
              ? 1.0
              : inner_radius * outer_radius / square(radius);
      return -extra_factor * point * one_over_distance_difference / radius;
    };

    for (auto falloff : {Falloff::Inverse, Falloff::Linear}) {
      expected.clear();
      expected.inner_sphericity = 1.0;
      expected.outer_sphericity = 1.0;
      expected.falloff = falloff;
      for (auto radius :
           {inner_radius, outer_radius, 0.5 * (inner_radius + outer_radius)}) {
        // Corners
        for (const double phi :
             {M_PI_4, 3.0 * M_PI_4, 5.0 * M_PI_4, 7.0 * M_PI_4}) {
          const std::array point =
              sph_to_cart(radius, acos(1.0 / sqrt(3.0)), phi);
          const double distortion =
              compute_distortion(make_not_null(&ylm), point, coefs);

          expected.points_to_check.emplace_back(point);
          const double transition =
              compute_transition(inner_radius, outer_radius, radius, falloff);
          expected.expected_transition.emplace_back(transition);
          expected.distortion.emplace_back(distortion);
          expected.expected_mapped_points.emplace_back(
              point * (1.0 - transition * distortion));
          expected.expected_gradient.emplace_back(
              compute_gradient(point, falloff));
        }

        // Along axes
        for (size_t i = 0; i < 3; i++) {
          for (double sign : {-1.0, 1.0}) {
            std::array<double, 3> point{0.0, 0.0, 0.0};
            gsl::at(point, i) = sign * radius;
            const double distortion =
                compute_distortion(make_not_null(&ylm), point, coefs);

            expected.points_to_check.emplace_back(point);
            const double transition =
                compute_transition(inner_radius, outer_radius, radius, falloff);
            expected.expected_transition.emplace_back(transition);
            expected.distortion.emplace_back(distortion);
            expected.expected_mapped_points.emplace_back(
                point * (1.0 - transition * distortion));
            expected.expected_gradient.emplace_back(
                compute_gradient(point, falloff));
          }
        }
      }

      test_points(expected, initial_time, functions_of_time);
    }
  }
  {
    // These are the corners of the cube
    INFO("Both cubes");
    // Gradient is less simple here
    const auto compute_gradient = [&](const std::array<double, 3>& point,
                                      const double inner_distance,
                                      const double outer_distance,
                                      const Falloff falloff) {
      const double radius = magnitude(point);
      const double one_over_denom = 1.0 / (outer_distance - inner_distance);
      size_t j = 0;
      double max = abs(point[0]);
      for (size_t i = 1; i < 3; i++) {
        if (abs(gsl::at(point, i)) > max) {
          j = i;
          max = abs(gsl::at(point, i));
        }
      }
      const size_t j_plus_1 = (j + 1) % 3;
      const size_t j_plus_2 = (j + 2) % 3;

      std::array<double, 3> grad_base =
          make_array<3, double>(1.0 / (sqrt(3.0) * radius * gsl::at(point, j)));
      gsl::at(grad_base, j) *=
          -(square(radius) - square(gsl::at(point, j))) / gsl::at(point, j);
      gsl::at(grad_base, j_plus_1) *= gsl::at(point, j_plus_1);
      gsl::at(grad_base, j_plus_2) *= gsl::at(point, j_plus_2);
      const std::array<double, 3> inner_grad = inner_radius * grad_base;
      const std::array<double, 3> outer_grad = outer_radius * grad_base;
      if (falloff == Falloff::Linear) {
        return (outer_grad - point / radius) * one_over_denom -
               (outer_distance - radius) * (outer_grad - inner_grad) *
                   square(one_over_denom);
      } else {
        const double a = -inner_distance * one_over_denom;
        const double b = -a * outer_distance;
        const std::array<double, 3> grad_a =
            (inner_distance * outer_grad - outer_distance * inner_grad) *
            square(one_over_denom);
        return grad_a - (a * outer_grad + outer_distance * grad_a) / radius -
               b * point / cube(radius);
      }
    };

    for (const auto falloff : {Falloff::Inverse, Falloff::Linear}) {
      CAPTURE(falloff);
      expected.clear();
      expected.inner_sphericity = 0.0;
      expected.outer_sphericity = 0.0;
      expected.falloff = falloff;

      // Corners
      for (const double radius :
           {inner_radius, outer_radius, 0.5 * (inner_radius + outer_radius)}) {
        for (const double phi :
             {M_PI_4, 3.0 * M_PI_4, 5.0 * M_PI_4, 7.0 * M_PI_4}) {
          const std::array point =
              sph_to_cart(radius, acos(1.0 / sqrt(3.0)), phi);
          const double distortion =
              compute_distortion(make_not_null(&ylm), point, coefs);

          expected.points_to_check.emplace_back(point);
          const double transition =
              compute_transition(inner_radius, outer_radius, radius, falloff);
          expected.expected_transition.emplace_back(transition);
          expected.distortion.emplace_back(distortion);
          expected.expected_mapped_points.emplace_back(
              point * (1.0 - transition * distortion));
          expected.expected_gradient.emplace_back(
              compute_gradient(point, inner_radius, outer_radius, falloff));
        }
      }

      // Along axes
      const double inner_distance = inner_radius / sqrt(3.0);
      const double outer_distance = outer_radius / sqrt(3.0);
      for (const double radius : {inner_distance, outer_distance,
                                  0.5 * (inner_distance + outer_distance)}) {
        for (size_t i = 0; i < 3; i++) {
          for (double sign : {-1.0, 1.0}) {
            std::array<double, 3> point{0.0, 0.0, 0.0};
            gsl::at(point, i) = sign * radius;
            const double distortion =
                compute_distortion(make_not_null(&ylm), point, coefs);

            expected.points_to_check.emplace_back(point);
            const double transition = compute_transition(
                inner_distance, outer_distance, radius, falloff);
            expected.expected_transition.emplace_back(transition);
            expected.distortion.emplace_back(distortion);
            expected.expected_mapped_points.emplace_back(
                point * (1.0 - transition * distortion));
            expected.expected_gradient.emplace_back(compute_gradient(
                point, inner_distance, outer_distance, falloff));
          }
        }
      }

      test_points(expected, initial_time, functions_of_time);
    }
  }
}

template <typename Generator>
void test_random_points(const gsl::not_null<Generator*> generator) {
  INFO("Test random points using shape map");
  std::uniform_real_distribution<double> coef_dist{-0.01, 0.01};
  std::uniform_int_distribution<size_t> num_dist{10_st, 20_st};
  const size_t num_points = num_dist(*generator);

  const double initial_time = 1.0;
  const size_t l_max = 4;
  const size_t num_coefs = 2 * square(l_max + 1);
  const std::array center{0.1, -0.2, 0.3};
  const std::string fot_name{"TheBean"};

  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      functions_of_time{};
  auto coefs = make_with_random_values<DataVector>(
      generator, make_not_null(&coef_dist), DataVector{num_coefs});
  functions_of_time[fot_name] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          initial_time, std::array{std::move(coefs)},
          std::numeric_limits<double>::infinity());

  // We'll keep the radii the same since they don't really matter
  const double inner_radius = 0.5;
  const double outer_radius = 10.0;

  std::uniform_real_distribution<double> angle_dist{0.0, 2.0 * M_PI};

  for (const auto& [inner_sphericity, outer_sphericity, falloff] :
       cartesian_product(std::array{1.0, 0.0}, std::array{1.0, 0.0},
                         std::array{Falloff::Inverse, Falloff::Linear})) {
    CAPTURE(inner_sphericity);
    CAPTURE(outer_sphericity);
    CAPTURE(falloff);
    // This guarantees the radius of the point is within the wedge
    std::uniform_real_distribution<double> radial_dist{
        inner_radius, outer_radius / sqrt(3.0)};

    std::unique_ptr<ShapeMapTransitionFunction> wedge =
        std::make_unique<Wedge>(inner_radius, outer_radius, inner_sphericity,
                                outer_sphericity, falloff);

    TimeDependent::Shape shape{center, l_max, l_max, std::move(wedge),
                               fot_name};

    for (size_t i = 0; i < num_points; i++) {
      const double radius = radial_dist(*generator);
      const double theta = 0.5 * angle_dist(*generator);
      const double phi = angle_dist(*generator);

      CAPTURE(radius);
      CAPTURE(theta);
      CAPTURE(phi);

      const std::array<double, 3> point =
          sph_to_cart(radius, theta, phi) + center;
      CAPTURE(point);

      const std::array<double, 3> mapped_point =
          shape(point, initial_time, functions_of_time);
      const std::optional<std::array<double, 3>> inverse_mapped_point =
          shape.inverse(mapped_point, initial_time, functions_of_time);

      CAPTURE(mapped_point);
      CAPTURE(inverse_mapped_point);

      CHECK(inverse_mapped_point.has_value());
      CHECK_ITERABLE_APPROX(point, inverse_mapped_point.value());
    }
  }
}

SPECTRE_TEST_CASE("Unit.Domain.CoordinateMaps.Shape.Wedge", "[Domain][Unit]") {
  MAKE_GENERATOR(generator);
  check_specific_points(make_not_null(&generator));
  test_random_points(make_not_null(&generator));
}
}  // namespace
}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
