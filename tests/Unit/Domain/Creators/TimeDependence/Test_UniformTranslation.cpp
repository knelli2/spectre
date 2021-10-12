// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <pup.h>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Domain/Creators/TimeDependence/TestHelpers.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace domain::creators::time_dependence {

namespace {
template <size_t MeshDim>
using Translation = domain::CoordinateMaps::TimeDependent::Translation<MeshDim>;

template <size_t MeshDim>
using ConcreteMap =
    domain::CoordinateMap<Frame::Grid, Frame::Inertial, Translation<MeshDim>>;

template <size_t MeshDim>
ConcreteMap<MeshDim> create_coord_map(const std::string& f_of_t_name) {
  return ConcreteMap<MeshDim>{Translation<MeshDim>{f_of_t_name}};
}

template <size_t MeshDim>
void test(const std::unique_ptr<TimeDependence<MeshDim>>& time_dep_unique_ptr,
          const double initial_time, const std::string& f_of_t_name) {
  MAKE_GENERATOR(gen);
  CAPTURE(initial_time);
  CAPTURE(f_of_t_name);

  CHECK_FALSE(time_dep_unique_ptr->is_none());

  // We downcast to the expected derived class to make sure that factory
  // creation worked correctly. In order to maximize code reuse this check is
  // done here as opposed to separately elsewhere.
  const auto* const time_dep = dynamic_cast<const UniformTranslation<MeshDim>*>(
      time_dep_unique_ptr.get());
  REQUIRE(time_dep != nullptr);

  // Test coordinate maps
  UniformCustomDistribution<size_t> dist_size_t{1, 10};
  const size_t num_blocks = dist_size_t(gen);
  CAPTURE(num_blocks);

  const auto expected_block_map = create_coord_map<MeshDim>(f_of_t_name);

  const auto block_maps = time_dep_unique_ptr->block_maps(num_blocks);
  for (const auto& block_map_unique_ptr : block_maps) {
    const auto* const block_map =
        dynamic_cast<const ConcreteMap<MeshDim>*>(block_map_unique_ptr.get());
    REQUIRE(block_map != nullptr);
    CHECK(*block_map == expected_block_map);
  }

  // Test functions of time
  const auto functions_of_time = time_dep_unique_ptr->functions_of_time();
  CHECK(functions_of_time.count(f_of_t_name) == 1);

  // Test map for composition
  CHECK(time_dep->map_for_composition() == expected_block_map);

  // For a random point at a random time check that the values agree. This is to
  // check that the internals were assigned the correct function of times.
  TIME_DEPENDENCE_GENERATE_COORDS(make_not_null(&gen), MeshDim, -1.0, 1.0);

  for (const auto& block_map : block_maps) {
    // We've checked equivalence above
    // (CHECK(*black_map == expected_block_map);), but have sometimes been
    // burned by incorrect operator== implementations so we check that the
    // mappings behave as expected.
    const double time_offset = dist(gen) + 1.2;
    CHECK_ITERABLE_APPROX(
        expected_block_map(grid_coords_dv, initial_time + time_offset,
                           functions_of_time),
        (*block_map)(grid_coords_dv, initial_time + time_offset,
                     functions_of_time));
    CHECK_ITERABLE_APPROX(
        expected_block_map(grid_coords_double, initial_time + time_offset,
                           functions_of_time),
        (*block_map)(grid_coords_double, initial_time + time_offset,
                     functions_of_time));

    CHECK_ITERABLE_APPROX(
        *expected_block_map.inverse(inertial_coords_double,
                                    initial_time + time_offset,
                                    functions_of_time),
        *block_map->inverse(inertial_coords_double, initial_time + time_offset,
                            functions_of_time));

    CHECK_ITERABLE_APPROX(
        expected_block_map.inv_jacobian(
            grid_coords_dv, initial_time + time_offset, functions_of_time),
        block_map->inv_jacobian(grid_coords_dv, initial_time + time_offset,
                                functions_of_time));
    CHECK_ITERABLE_APPROX(
        expected_block_map.inv_jacobian(
            grid_coords_double, initial_time + time_offset, functions_of_time),
        block_map->inv_jacobian(grid_coords_double, initial_time + time_offset,
                                functions_of_time));

    CHECK_ITERABLE_APPROX(
        expected_block_map.jacobian(grid_coords_dv, initial_time + time_offset,
                                    functions_of_time),
        block_map->jacobian(grid_coords_dv, initial_time + time_offset,
                            functions_of_time));
    CHECK_ITERABLE_APPROX(
        expected_block_map.jacobian(
            grid_coords_double, initial_time + time_offset, functions_of_time),
        block_map->jacobian(grid_coords_double, initial_time + time_offset,
                            functions_of_time));
  }
}

void test_equivalence() {
  {
    UniformTranslation<1> ut0{1.0, {{2.0}}};
    UniformTranslation<1> ut1{1.2, {{2.0}}};
    UniformTranslation<1> ut2{1.0, {{3.0}}};
    CHECK(ut0 == ut0);
    CHECK_FALSE(ut0 != ut0);
    CHECK(ut0 != ut1);
    CHECK_FALSE(ut0 == ut1);
    CHECK(ut0 != ut2);
    CHECK_FALSE(ut0 == ut2);
  }
  {
    UniformTranslation<2> ut0{1.0, {{2.0, 4.0}}};
    UniformTranslation<2> ut1{1.2, {{2.0, 4.0}}};
    UniformTranslation<2> ut2{1.0, {{3.0, 4.0}}};
    UniformTranslation<2> ut3{1.0, {{2.0, 5.0}}};
    CHECK(ut0 == ut0);
    CHECK_FALSE(ut0 != ut0);
    CHECK(ut0 != ut1);
    CHECK_FALSE(ut0 == ut1);
    CHECK(ut0 != ut2);
    CHECK_FALSE(ut0 == ut2);
    CHECK(ut0 != ut3);
    CHECK_FALSE(ut0 == ut3);
  }
  {
    UniformTranslation<3> ut0{1.0, {{2.0, 4.0, 6.0}}};
    UniformTranslation<3> ut1{1.2, {{2.0, 4.0, 6.0}}};
    UniformTranslation<3> ut2{1.0, {{3.0, 4.0, 6.0}}};
    UniformTranslation<3> ut3{1.0, {{2.0, 5.0, 6.0}}};
    UniformTranslation<3> ut4{1.0, {{2.0, 4.0, 7.0}}};
    CHECK(ut0 == ut0);
    CHECK_FALSE(ut0 != ut0);
    CHECK(ut0 != ut1);
    CHECK_FALSE(ut0 == ut1);
    CHECK(ut0 != ut2);
    CHECK_FALSE(ut0 == ut2);
    CHECK(ut0 != ut3);
    CHECK_FALSE(ut0 == ut3);
    CHECK(ut0 != ut4);
    CHECK_FALSE(ut0 == ut4);
  }
}

SPECTRE_TEST_CASE("Unit.Domain.Creators.TimeDependence.UniformTranslation",
                  "[Domain][Unit]") {
  const double initial_time = 1.3;

  {
    // 1d
    const std::array<double, 1> velocity{{2.4}};
    const std::string f_of_t_name =
        "UniformTranslation::vel=" + get_output(velocity) +
        "::t_0=" + get_output(initial_time);
    const std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
        time_dep =
            std::make_unique<UniformTranslation<1>>(initial_time, velocity);
    test(time_dep, initial_time, f_of_t_name);
    test(time_dep->get_clone(), initial_time, f_of_t_name);

    test(TestHelpers::test_creation<std::unique_ptr<TimeDependence<1>>>(
             "UniformTranslation:\n"
             "  InitialTime: 1.3\n"
             "  Velocity: [2.4]\n"),
         initial_time, f_of_t_name);
  }

  {
    // 2d
    const std::array<double, 2> velocity{{2.4, 3.1}};
    const std::string f_of_t_name =
        "UniformTranslation::vel=" + get_output(velocity) +
        "::t_0=" + get_output(initial_time);
    const std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
        time_dep =
            std::make_unique<UniformTranslation<2>>(initial_time, velocity);
    test(time_dep, initial_time, f_of_t_name);
    test(time_dep->get_clone(), initial_time, f_of_t_name);

    test(TestHelpers::test_creation<std::unique_ptr<TimeDependence<2>>>(
             "UniformTranslation:\n"
             "  InitialTime: 1.3\n"
             "  Velocity: [2.4, 3.1]\n"),
         initial_time, f_of_t_name);
  }

  {
    // 3d
    const std::array<double, 3> velocity{{2.4, 3.1, -1.2}};
    const std::string f_of_t_name =
        "UniformTranslation::vel=" + get_output(velocity) +
        "::t_0=" + get_output(initial_time);
    const std::unique_ptr<domain::creators::time_dependence::TimeDependence<3>>
        time_dep =
            std::make_unique<UniformTranslation<3>>(initial_time, velocity);
    test(time_dep, initial_time, f_of_t_name);
    test(time_dep->get_clone(), initial_time, f_of_t_name);

    test(TestHelpers::test_creation<std::unique_ptr<TimeDependence<3>>>(
             "UniformTranslation:\n"
             "  InitialTime: 1.3\n"
             "  Velocity: [2.4, 3.1, -1.2]\n"),
         initial_time, f_of_t_name);
  }

  test_equivalence();
}
}  // namespace

}  // namespace domain::creators::time_dependence
