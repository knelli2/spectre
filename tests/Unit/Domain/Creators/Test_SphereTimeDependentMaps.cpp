// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>

#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/SphereTimeDependentMaps.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Framework/TestCreation.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"

namespace Frame {
struct Grid;
struct Distorted;
struct Inertial;
}  // namespace Frame

namespace domain::creators::sphere {
namespace {
std::string create_option_string(const bool use_non_zero_shape) {
  return "InitialTime: 1.5\n"
         "ShapeMap:\n"
         "  LMax: 8\n" +
         (use_non_zero_shape ? "  InitialValues:\n"
                               "    Mass: 1.0\n"
                               "    Spin: [0.0, 0.0, 0.0]\n"s
                             : "  InitialValues: Spherical\n"s) +
         "RotationMap: None\n"
         "ExpansionMap: None\n"
         "TranslationMap:\n"
         "  InitialValues: [[0.1, -3.2, 1.1], [-0.3, 0.5, -0.7],"
         " [0.1, -0.4, 0.02]]\n"s;
}
void test(const bool use_non_zero_shape) {
  CAPTURE(use_non_zero_shape);
  const double initial_time = 1.5;
  const size_t l_max = 8;

  auto time_dep_options = TestHelpers::test_creation<TimeDependentMapOptions>(
      create_option_string(use_non_zero_shape));

  std::unordered_map<std::string, double> expiration_times{
      {TimeDependentMapOptions::shape_name, 15.5},
      {TimeDependentMapOptions::size_name,
       std::numeric_limits<double>::infinity()},
      {TimeDependentMapOptions::translation_name,
       std::numeric_limits<double>::infinity()}};

  // These are hard-coded so this is just a regression test
  CHECK(TimeDependentMapOptions::size_name == "Size"s);
  CHECK(TimeDependentMapOptions::shape_name == "Shape"s);
  CHECK(TimeDependentMapOptions::rotation_name == "Rotation"s);
  CHECK(TimeDependentMapOptions::expansion_name == "Expansion"s);
  CHECK(TimeDependentMapOptions::expansion_outer_boundary_name ==
        "ExpansionOuterBoundary"s);
  CHECK(TimeDependentMapOptions::translation_name == "Translation"s);

  using PP2 = domain::FunctionsOfTime::PiecewisePolynomial<2>;
  using PP3 = domain::FunctionsOfTime::PiecewisePolynomial<3>;
  DataVector size_func{1, 0.0};
  DataVector size_deriv{1, 0.0};
  DataVector size_2nd_deriv{1, 0.0};
  PP3 size{
      initial_time,
      std::array<DataVector, 4>{{size_func, size_deriv, size_2nd_deriv, {0.0}}},
      expiration_times.at(TimeDependentMapOptions::size_name)};
  const DataVector shape_zeros{ylm::Spherepack::spectral_size(l_max, l_max),
                               0.0};
  PP2 shape_all_zero{
      initial_time,
      std::array<DataVector, 3>{shape_zeros, shape_zeros, shape_zeros},
      expiration_times.at(TimeDependentMapOptions::shape_name)};

  DataVector initial_translation_center{{0.1, -3.2, 1.1}};
  DataVector initial_velocity{{-0.3, 0.5, -0.7}};
  DataVector initial_translation_acceleration{{0.1, -0.4, 0.02}};
  PP2 translation_non_zero{
      initial_time,
      std::array<DataVector, 3>{{initial_translation_center, initial_velocity,
                                 initial_translation_acceleration}},
      expiration_times.at(TimeDependentMapOptions::translation_name)};

  const std::array<double, 3> center{-5.0, -0.01, -0.02};
  const double inner_radius = 0.5;
  const double outer_radius = 2.1;
  const double transition_inner_radius = 3.1;
  const double transition_outer_radius = 3.9;

  const auto functions_of_time =
      time_dep_options.create_functions_of_time(inner_radius, expiration_times);

  CHECK(dynamic_cast<PP3&>(
            *functions_of_time.at(TimeDependentMapOptions::size_name).get()) ==
        size);
  // This will always be zero based on our options above. Just easier to check
  // that way
  CHECK(dynamic_cast<PP2&>(
            *functions_of_time.at(TimeDependentMapOptions::shape_name).get()) ==
        shape_all_zero);

  CHECK(dynamic_cast<PP2&>(
            *functions_of_time.at(TimeDependentMapOptions::translation_name)
                 .get()) == translation_non_zero);

  for (const bool include_distorted : make_array(true, false)) {
    for (const bool use_rigid : make_array(true, false)) {
      time_dep_options.build_maps(
          center, std::pair<double, double>{inner_radius, outer_radius},
          std::pair<double, double>{transition_inner_radius,
                                    transition_outer_radius});

      const auto grid_to_distorted_map =
          time_dep_options.grid_to_distorted_map(include_distorted);
      const auto grid_to_inertial_map =
          time_dep_options.grid_to_inertial_map(include_distorted, use_rigid);
      const auto distorted_to_inertial_map =
          time_dep_options.distorted_to_inertial_map(include_distorted);

      // All of these maps are tested individually. Rather than going through
      // the effort of coming up with a source coordinate and calculating
      // analytically what we would get after it's mapped, we just check that
      // whether it's supposed to be a nullptr and if it's not, if it's the
      // identity and that the jacobians are time dependent.
      const auto check_map = [](const auto& map, const bool is_null,
                                const bool is_identity) {
        if (is_null) {
          CHECK(map == nullptr);
        } else {
          CHECK(map->is_identity() == is_identity);
          CHECK(map->inv_jacobian_is_time_dependent() == not is_identity);
          CHECK(map->jacobian_is_time_dependent() == not is_identity);
        }
      };

      // There is no null pointer in the grid to inertial map
      check_map(grid_to_inertial_map, false, false);

      check_map(grid_to_distorted_map, not include_distorted, false);

      // If no shape distortion, there is only the rotation, expansion and
      // translation maps
      check_map(distorted_to_inertial_map, not include_distorted, false);
    }
  }
}

void test_shape_initial_values() {
  KerrSchildFromBoyerLindquist shape_params{1.4, {0.1, 0.2, -0.3}};

  auto option_shape_params =
      TestHelpers::test_creation<KerrSchildFromBoyerLindquist>(
          "Mass: 1.4\n"
          "Spin: [0.1, 0.2, -0.3]");

  CHECK(shape_params.mass == option_shape_params.mass);
  CHECK(shape_params.spin == option_shape_params.spin);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.Creators.SphereTimeDependentMaps",
                  "[Domain][Unit]") {
  test_shape_initial_values();

  test(true);
  test(false);
}
}  // namespace domain::creators::sphere
