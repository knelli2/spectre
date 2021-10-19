// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <memory>
#include <random>
#include <vector>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.tpp"
#include "Domain/Creators/TimeDependence/CompositionCubicScaleAndUniformRotationAboutZAxis.hpp"
#include "Domain/Creators/TimeDependence/CubicScale.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformRotationAboutZAxis.hpp"
#include "Domain/FunctionsOfTime/CombineFunctionsOfTime.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Domain/Creators/TimeDependence/TestHelpers.hpp"

namespace domain::creators::time_dependence {
namespace {
template <size_t Dim>
using CubicScaleMapTimeDep =
    domain::CoordinateMaps::TimeDependent::CubicScale<Dim>;
using RotationMapTimeDep = domain::CoordinateMaps::TimeDependent::Rotation<2>;
using Identity = domain::CoordinateMaps::Identity<1>;
template <size_t Dim>
using RotationAboutZAxis =
    tmpl::conditional_t<Dim == 2, RotationMapTimeDep,
                        domain::CoordinateMaps::TimeDependent::ProductOf2Maps<
                            RotationMapTimeDep, Identity>>;

template <size_t Dim>
using CubicScaleMap = domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                            CubicScaleMapTimeDep<Dim>>;

template <size_t Dim>
using RotationMap = domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                          RotationAboutZAxis<Dim>>;

template <size_t Dim>
using CompositionMap =
    domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                          CubicScaleMapTimeDep<Dim>, RotationAboutZAxis<Dim>>;

template <typename T0, size_t Dim>
void test_impl(
    const gsl::not_null<T0> gen, const double initial_time,
    const std::unique_ptr<TimeDependence<Dim>>& time_dep,
    const std::vector<std::unique_ptr<TimeDependence<Dim>>>& expected_time_deps,
    const std::vector<std::string>& expected_f_of_t_names) {
  // Test coordinate maps
  UniformCustomDistribution<size_t> dist_size_t{1, 10};
  const size_t num_blocks = dist_size_t(*gen);
  CAPTURE(num_blocks);

  CHECK_FALSE(time_dep->is_none());

  // Test functions of time without expiration times
  {
    const auto functions_of_time = time_dep->functions_of_time();
    auto functions_of_time_for_expected =
        expected_time_deps[0]->functions_of_time();
    for (size_t i = 1; i < expected_time_deps.size(); i++) {
      functions_of_time_for_expected =
          domain::FunctionsOfTime::combine_functions_of_time(
              functions_of_time_for_expected,
              expected_time_deps[i]->functions_of_time());
    }

    REQUIRE(functions_of_time.size() == expected_f_of_t_names.size());
    for (const auto& f_of_t_name : expected_f_of_t_names) {
      CHECK(functions_of_time.count(f_of_t_name) == 1);
      CHECK(functions_of_time.at(f_of_t_name)->time_bounds()[1] ==
            std::numeric_limits<double>::infinity());
    }
  }
  // Test functions of time with expiration times
  {
    const double init_expr_time = 5.0;
    std::unordered_map<std::string, double> init_expr_times{};
    for (auto& name : expected_f_of_t_names) {
      init_expr_times[name] = init_expr_time;
    }
    const auto functions_of_time = time_dep->functions_of_time(init_expr_times);
    auto functions_of_time_for_expected =
        expected_time_deps[0]->functions_of_time(init_expr_times);
    for (size_t i = 1; i < expected_time_deps.size(); i++) {
      functions_of_time_for_expected =
          domain::FunctionsOfTime::combine_functions_of_time(
              functions_of_time_for_expected,
              expected_time_deps[i]->functions_of_time(init_expr_times));
    }

    REQUIRE(functions_of_time.size() == expected_f_of_t_names.size());
    for (const auto& f_of_t_name : expected_f_of_t_names) {
      CHECK(functions_of_time.count(f_of_t_name) == 1);
      CHECK(functions_of_time.at(f_of_t_name)->time_bounds()[1] ==
            init_expr_time);
    }
  }

  const auto functions_of_time = time_dep->functions_of_time();
  auto functions_of_time_for_expected =
      expected_time_deps[0]->functions_of_time();
  for (size_t i = 1; i < expected_time_deps.size(); i++) {
    functions_of_time_for_expected =
        domain::FunctionsOfTime::combine_functions_of_time(
            functions_of_time_for_expected,
            expected_time_deps[i]->functions_of_time());
  }

  // For a random point at a random time check that the values agree. This is to
  // check that the internals were assigned the correct function of times.
  TIME_DEPENDENCE_GENERATE_COORDS(gen, Dim, -1.0, 1.0);

  const auto block_maps = time_dep->block_maps(num_blocks);
  auto expected_block_maps = expected_time_deps[0]->block_maps(num_blocks);
  for (size_t i = 1; i < expected_time_deps.size(); i++) {
    auto temp_block_maps = expected_time_deps[i]->block_maps(num_blocks);
    for (size_t j = 0; j < temp_block_maps.size(); j++) {
      expected_block_maps[j] =
          std::make_unique<CompositionMap<Dim>>(domain::push_back(
              dynamic_cast<CubicScaleMap<Dim>&>(*(expected_block_maps[j])),
              dynamic_cast<RotationMap<Dim>&>(*(temp_block_maps[j]))));
    }
  }

  REQUIRE(block_maps.size() == expected_block_maps.size());

  for (size_t i = 0; i < block_maps.size(); ++i) {
    const double time_offset = dist(*gen) + 1.2;
    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]->operator()(grid_coords_dv,
                                           initial_time + time_offset,
                                           functions_of_time_for_expected),
        block_maps[i]->operator()(grid_coords_dv, initial_time + time_offset,
                                  functions_of_time));
    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]->operator()(grid_coords_double,
                                           initial_time + time_offset,
                                           functions_of_time_for_expected),
        block_maps[i]->operator()(
            grid_coords_double, initial_time + time_offset, functions_of_time));

    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]
            ->inverse(inertial_coords_double, initial_time + time_offset,
                      functions_of_time_for_expected)
            .value(),
        block_maps[i]
            ->inverse(inertial_coords_double, initial_time + time_offset,
                      functions_of_time)
            .value());

    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]->inv_jacobian(grid_coords_dv,
                                             initial_time + time_offset,
                                             functions_of_time_for_expected),
        block_maps[i]->inv_jacobian(grid_coords_dv, initial_time + time_offset,
                                    functions_of_time));
    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]->inv_jacobian(grid_coords_double,
                                             initial_time + time_offset,
                                             functions_of_time_for_expected),
        block_maps[i]->inv_jacobian(
            grid_coords_double, initial_time + time_offset, functions_of_time));

    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]->jacobian(grid_coords_dv,
                                         initial_time + time_offset,
                                         functions_of_time_for_expected),
        block_maps[i]->jacobian(grid_coords_dv, initial_time + time_offset,
                                functions_of_time));
    CHECK_ITERABLE_APPROX(
        expected_block_maps[i]->jacobian(grid_coords_double,
                                         initial_time + time_offset,
                                         functions_of_time_for_expected),
        block_maps[i]->jacobian(grid_coords_double, initial_time + time_offset,
                                functions_of_time));
  }
}

template <size_t Dim, typename T>
void test_composition_cubic_scale_uniform_rotation(
    const gsl::not_null<T> gen, const double initial_time,
    const bool use_linear_scaling) {
  static_assert(Dim != 1,
                "CompositionCubicScaleAndUniformRotationAboutZAxis undefined "
                "for Dim == 1.");
  INFO("Test composition of cubic scale and uniform rotation");
  using Composition = CompositionCubicScaleAndUniformRotationAboutZAxis<Dim>;
  UniformCustomDistribution<double> dist_double{0.0, 1.0};
  const double outer_boundary = 10.4;
  const std::array<double, 2> initial_expansion{{1.0, 1.0}};
  const std::array<double, 2> velocity{{-0.1, 0.0}};
  const std::array<double, 2> acceleration{{-0.05, 0.0}};
  const std::string cs_name0{"CubicScale"s + (use_linear_scaling ? "" : "A")};
  const std::string cs_name1{"CubicScale"s + (use_linear_scaling ? "" : "B")};
  std::array<std::string, 2> f_of_t_names{cs_name0, cs_name1};

  CubicScale<Dim> time_dep0(initial_time, outer_boundary, use_linear_scaling,
                            initial_expansion, velocity, acceleration);

  const double angular_velocity = 2.4;
  const std::string rot_name{"Rotation"};

  UniformRotationAboutZAxis<Dim> time_dep1(initial_time, angular_velocity);

  std::vector<std::unique_ptr<TimeDependence<Dim>>> expected_time_deps{};
  expected_time_deps.push_back(time_dep0.get_clone());
  expected_time_deps.push_back(time_dep1.get_clone());

  std::vector<std::string> expected_names{rot_name, cs_name0};
  if (not use_linear_scaling) {
    expected_names.push_back(cs_name1);
  }

  const std::unique_ptr<TimeDependence<Dim>> time_dep =
      std::make_unique<Composition>(time_dep0.get_clone(),
                                    time_dep1.get_clone());

  test_impl(gen, initial_time, time_dep, expected_time_deps, expected_names);

  test_impl(gen, initial_time, time_dep->get_clone(), expected_time_deps,
            expected_names);
}

template <typename T>
void test_options(const gsl::not_null<T> gen, const double initial_time,
                  const bool use_linear_scaling) {
  INFO("Test create by options");
  const double outer_boundary = 10.4;
  const std::array<double, 2> initial_expansion{{1.0, 1.0}};
  const std::array<double, 2> velocity{{-0.1, 0.0}};
  const std::array<double, 2> acceleration{{-0.05, 0.0}};
  const std::string cs_name0{"CubicScale"s + (use_linear_scaling ? "" : "A")};
  const std::string cs_name1{"CubicScale"s + (use_linear_scaling ? "" : "B")};
  std::array<std::string, 2> f_of_t_names{cs_name0, cs_name1};

  CubicScale<2> time_dep0(initial_time, outer_boundary, use_linear_scaling,
                          initial_expansion, velocity, acceleration);

  const double angular_velocity = 2.4;
  const std::string rot_name{"Rotation"};

  UniformRotationAboutZAxis<2> time_dep1(initial_time, angular_velocity);

  std::vector<std::unique_ptr<TimeDependence<2>>> expected_time_deps{};
  expected_time_deps.push_back(time_dep0.get_clone());
  expected_time_deps.push_back(time_dep1.get_clone());

  std::vector<std::string> expected_names{rot_name, cs_name0};
  if (not use_linear_scaling) {
    expected_names.push_back(cs_name1);
  }

  const std::string linear_scale_str =
      "      UseLinearScaling: "s + (use_linear_scaling ? "true\n" : "false\n");

  const auto created_with_options =
      TestHelpers::test_creation<std::unique_ptr<TimeDependence<2>>>(
          "CompositionCubicScaleAndUniformRotationAboutZAxis:\n"
          "  CubicScale:\n"
          "    CubicScale:\n"
          "      InitialTime: 1.3\n"
          "      OuterBoundary: 10.4\n"
          "      InitialExpansion: [1.0, 1.0]\n" +
          linear_scale_str +
          "      Velocity: [-0.1, 0.0]\n"
          "      Acceleration: [-0.05, 0.0]\n"
          "  UniformRotationAboutZAxis:\n"
          "    UniformRotationAboutZAxis:\n"
          "      InitialTime: 1.3\n"
          "      AngularVelocity: 2.4\n");

  test_impl(gen, initial_time, created_with_options, expected_time_deps,
            expected_names);
}

SPECTRE_TEST_CASE(
    "Unit.Domain.Creators.TimeDependence.CompositionCubicScaleUniformRotation",
    "[Domain][Unit]") {
  MAKE_GENERATOR(gen);
  const double initial_time = 1.3;
  test_composition_cubic_scale_uniform_rotation<2>(make_not_null(&gen),
                                                   initial_time, true);
  test_composition_cubic_scale_uniform_rotation<2>(make_not_null(&gen),
                                                   initial_time, false);
  test_composition_cubic_scale_uniform_rotation<3>(make_not_null(&gen),
                                                   initial_time, true);
  test_composition_cubic_scale_uniform_rotation<3>(make_not_null(&gen),
                                                   initial_time, false);
  test_options(make_not_null(&gen), initial_time, true);
  test_options(make_not_null(&gen), initial_time, false);
}
}  // namespace
}  // namespace domain::creators::time_dependence
