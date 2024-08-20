// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <pup.h>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/CoordinateMaps/Rotation.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/CoordinateMaps/Wedge.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Domain/Tags.hpp"
#include "Framework/TestHelpers.hpp"

namespace domain {
namespace {
using DV = DataVector;
using Affine = CoordinateMaps::Affine;

template <size_t Dim, typename S, typename T, typename U>
void test_element_impl(
    bool test_inverse, const ElementId<Dim>& element_id, const S& affine_map,
    const T& first_map, const U& second_map,
    const tnsr::I<double, Dim, Frame::ElementLogical>& logical_point_double,
    const tnsr::I<DV, Dim, Frame::ElementLogical>& logical_point_dv) {
  PUPable_reg(
      SINGLE_ARG(CoordinateMap<Frame::BlockLogical, Frame::Inertial, T, U>));
  const auto composed_map =
      make_coordinate_map<Frame::ElementLogical, Frame::Inertial>(
          affine_map, first_map, second_map);

  ElementMap element_map{
      element_id,
      make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
          first_map, second_map)};
  ElementMap element_map_deserialized = serialize_and_deserialize(element_map);

  CHECK(element_map(logical_point_dv) == composed_map(logical_point_dv));
  CHECK(element_map(logical_point_double) ==
        composed_map(logical_point_double));
  CHECK(element_map_deserialized(logical_point_dv) ==
        composed_map(logical_point_dv));
  CHECK(element_map_deserialized(logical_point_double) ==
        composed_map(logical_point_double));

  const tnsr::I<double, Dim, Frame::Inertial> inertial_point_double =
      composed_map(logical_point_double);
  const tnsr::I<DV, Dim, Frame::Inertial> inertial_point_dv =
      composed_map(logical_point_dv);

  if (test_inverse) {
    CHECK(element_map.inverse(inertial_point_double) ==
          composed_map.inverse(inertial_point_double).value());
    CHECK(element_map_deserialized.inverse(inertial_point_double) ==
          composed_map.inverse(inertial_point_double).value());
  }

  CHECK_ITERABLE_APPROX(element_map.inv_jacobian(logical_point_dv),
                        composed_map.inv_jacobian(logical_point_dv));
  CHECK_ITERABLE_APPROX(element_map.inv_jacobian(logical_point_double),
                        composed_map.inv_jacobian(logical_point_double));
  CHECK_ITERABLE_APPROX(element_map_deserialized.inv_jacobian(logical_point_dv),
                        composed_map.inv_jacobian(logical_point_dv));
  CHECK_ITERABLE_APPROX(
      element_map_deserialized.inv_jacobian(logical_point_double),
      composed_map.inv_jacobian(logical_point_double));

  CHECK_ITERABLE_APPROX(element_map.inv_jacobian(logical_point_dv),
                        composed_map.inv_jacobian(logical_point_dv));
  CHECK_ITERABLE_APPROX(element_map.jacobian(logical_point_double),
                        composed_map.jacobian(logical_point_double));
  CHECK_ITERABLE_APPROX(element_map_deserialized.inv_jacobian(logical_point_dv),
                        composed_map.inv_jacobian(logical_point_dv));
  CHECK_ITERABLE_APPROX(element_map_deserialized.jacobian(logical_point_double),
                        composed_map.jacobian(logical_point_double));

  CHECK(element_map.block_map() ==
        *(make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
            first_map, second_map)));
}

template <size_t Dim>
void test_element_map();

template <>
void test_element_map<1>() {
  auto segment_ids = std::array<SegmentId, 1>({{SegmentId(2, 3)}});
  ElementId<1> element_id(0, segment_ids);
  const tnsr::I<double, 1, Frame::ElementLogical> logical_point_double{0.0};
  const tnsr::I<DV, 1, Frame::ElementLogical> logical_point_dv(
      DV{-1.0, -0.5, 0.0, 0.5, 1.0});

  const Affine affine_map{-1.0, 1.0, 0.5, 1.0};
  const Affine first_map{-1.0, 1.0, 2.0, 8.0};

  // Aligned maps
  test_element_impl(true, element_id, affine_map, first_map,
                    Affine{2.0, 8.0, -2.0, -1.0}, logical_point_double,
                    logical_point_dv);

  // Flip axis in second map
  test_element_impl(true, element_id, affine_map, first_map,
                    Affine{2.0, 8.0, 2.0, -1.0}, logical_point_double,
                    logical_point_dv);
}

template <>
void test_element_map<2>() {
  using Rotate = CoordinateMaps::Rotation<2>;
  using Wedge2D = CoordinateMaps::Wedge<2>;

  auto segment_ids =
      std::array<SegmentId, 2>({{SegmentId(2, 3), SegmentId(1, 0)}});
  ElementId<2> element_id(0, segment_ids);

  const tnsr::I<double, 2, Frame::ElementLogical> logical_point_double(
      std::array<double, 2>{{-0.5, 0.5}});
  const tnsr::I<DV, 2, Frame::ElementLogical> logical_point_dv(
      std::array<DV, 2>{
          {DV{-1.0, 0.5, 0.0, 0.5, 1.0}, DV{-1.0, 0.5, 0.0, 0.5, 1.0}}});

  const CoordinateMaps::ProductOf2Maps<Affine, Affine> affine_map(
      Affine{-1.0, 1.0, 0.5, 1.0}, Affine{-1.0, 1.0, -1.0, 0.0});
  const auto first_map = Rotate(2.);

  // Test with two rotations
  test_element_impl(true, element_id, affine_map, first_map, Rotate(1.8472),
                    logical_point_double, logical_point_dv);

  // test with a rotation and a wedge
  test_element_impl(
      false, element_id, affine_map, first_map,
      Wedge2D(3., 7., 0.0, 1.0, 1.0, {{0.0, 0.0}},
              OrientationMap<2>{std::array<Direction<2>, 2>{
                  {Direction<2>::lower_xi(), Direction<2>::lower_eta()}}},
              false),
      logical_point_double, logical_point_dv);
}

template <>
void test_element_map<3>() {
  using Rotate = CoordinateMaps::Rotation<3>;

  auto segment_ids = std::array<SegmentId, 3>(
      {{SegmentId(2, 3), SegmentId(1, 0), SegmentId(2, 1)}});
  ElementId<3> element_id(0, segment_ids);

  const tnsr::I<double, 3, Frame::ElementLogical> logical_point_double{
      {{0.78, 0.12, -0.37}}};
  const tnsr::I<DataVector, 3, Frame::ElementLogical> logical_point_dv{
      {{DV{-1.0, 0.5, 0.0, 0.5, 1.0}, DV{-1.0, 0.5, 0.0, 0.5, 1.0},
        DV{-1.0, 0.5, 0.0, 0.5, 1.0}}}};

  const CoordinateMaps::ProductOf3Maps<Affine, Affine, Affine> affine_map(
      Affine{-1.0, 1.0, 0.5, 1.0}, Affine{-1.0, 1.0, -1.0, 0.0},
      Affine{-1.0, 1.0, -0.5, 0.0});
  const auto first_map = Rotate{M_PI_4, M_PI_4, M_PI_2};

  // test with 2 rotations
  test_element_impl(true, element_id, affine_map, first_map,
                    Rotate{M_PI_2, M_PI_4, M_PI_4}, logical_point_double,
                    logical_point_dv);

  // test with rotation and wedge
  test_element_impl(
      true, element_id, affine_map, first_map,
      CoordinateMaps::Wedge<3>{3.0,
                               7.0,
                               0.8,
                               0.9,
                               1.0,
                               {{0.0, 0.0, 0.0}},
                               OrientationMap<3>::create_aligned(),
                               true},
      logical_point_double, logical_point_dv);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.ElementMap", "[Unit][Domain]") {
  test_element_map<1>();
  test_element_map<2>();
  test_element_map<3>();

  {
    INFO("Test with Block");
    const ElementId<1> element_id{0, {{{2, 3}}}};
    const Affine block_map{-1.0, 1.0, 2.0, 8.0};
    const tnsr::I<double, 1, Frame::ElementLogical> xi{{{1.}}};
    Block<1> block{
        make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
            block_map),
        0,
        {}};
    {
      INFO("Time-independent");
      {
        INFO("Grid frame");
        ElementMap<1, Frame::Grid> element_map{element_id, block};
        ElementMap<1, Frame::Grid> expected_element_map{
            element_id,
            make_coordinate_map_base<Frame::BlockLogical, Frame::Grid>(
                block_map)};
        CHECK(element_map(xi) == expected_element_map(xi));
      }
      {
        INFO("Inertial frame");
        ElementMap<1, Frame::Inertial> element_map{element_id, block};
        ElementMap<1, Frame::Inertial> expected_element_map{
            element_id,
            make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
                block_map)};
        CHECK(element_map(xi) == expected_element_map(xi));
      }
    }
    {
      INFO("Time-dependent");
      const domain::CoordinateMaps::TimeDependent::Translation<1>
          grid_to_inertial_map{"Translation"};
      block.inject_time_dependent_map(
          make_coordinate_map_base<Frame::Grid, Frame::Inertial>(
              grid_to_inertial_map));
      {
        INFO("Grid frame");
        ElementMap<1, Frame::Grid> element_map{element_id, block};
        ElementMap<1, Frame::Grid> expected_element_map{
            element_id,
            make_coordinate_map_base<Frame::BlockLogical, Frame::Grid>(
                block_map)};
        CHECK(element_map(xi) == expected_element_map(xi));
      }
      {
        INFO("Inertial frame");
        ElementMap<1, Frame::Inertial> element_map{element_id, block};
        ElementMap<1, Frame::Inertial> expected_element_map{
            element_id,
            make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
                block_map, grid_to_inertial_map)};
        const double time = 1.;
        domain::FunctionsOfTimeMap functions_of_time{};
        functions_of_time["Translation"] =
            std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<1>>(
                time, std::array<DataVector, 2>{{{1, 2.}, {1, 0.}}},
                std::numeric_limits<double>::infinity());
        CHECK(element_map(xi, time, functions_of_time) ==
              expected_element_map(xi, time, functions_of_time));
        CHECK(get<0>(element_map(xi, time, functions_of_time)) == 10.);
      }
    }
  }
}
}  // namespace domain
