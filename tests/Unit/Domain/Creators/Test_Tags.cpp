// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>
#include <type_traits>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/BinaryCompactObject.hpp"
#include "Domain/Creators/Brick.hpp"
#include "Domain/Creators/OptionTags.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/ExternalBoundaryConditions.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/Creators/Tags/InitialExtents.hpp"
#include "Domain/Creators/Tags/InitialRefinementLevels.hpp"
#include "Domain/Creators/Tags/ObjectCenter.hpp"
#include "Domain/FunctionsOfTime/OptionTags.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"

#include <memory>
#include <unordered_map>

#include "Domain/CoordinateMaps/TimeDependent/Shape.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/CoordinateMaps/TimeDependent/SphericalCompression.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/YlmSpherepack.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain {
namespace {
template <size_t Dim>
void test_simple_tags() {
  TestHelpers::db::test_simple_tag<Tags::Domain<Dim>>("Domain");
  TestHelpers::db::test_simple_tag<Tags::InitialExtents<Dim>>("InitialExtents");
  TestHelpers::db::test_simple_tag<Tags::InitialRefinementLevels<Dim>>(
      "InitialRefinementLevels");
  TestHelpers::db::test_simple_tag<Tags::ExternalBoundaryConditions<Dim>>(
      "ExternalBoundaryConditions");
}

void test_center_tags() {
  TestHelpers::db::test_base_tag<Tags::ObjectCenter<ObjectLabel::A>>(
      "ObjectCenter");
  TestHelpers::db::test_base_tag<Tags::ObjectCenter<ObjectLabel::B>>(
      "ObjectCenter");
  TestHelpers::db::test_simple_tag<Tags::ExcisionCenter<ObjectLabel::A>>(
      "CenterObjectA");
  TestHelpers::db::test_simple_tag<Tags::ExcisionCenter<ObjectLabel::B>>(
      "CenterObjectB");

  using Object = domain::creators::BinaryCompactObject::Object;

  const std::unique_ptr<DomainCreator<3>> domain_creator =
      std::make_unique<domain::creators::BinaryCompactObject>(
          Object{0.2, 5.0, 8.0, true, true}, Object{0.6, 4.0, -5.5, true, true},
          100.0, 500.0, 1_st, 5_st);

  const auto grid_center_A =
      Tags::ExcisionCenter<ObjectLabel::A>::create_from_options(domain_creator);
  const auto grid_center_B =
      Tags::ExcisionCenter<ObjectLabel::B>::create_from_options(domain_creator);

  CHECK(grid_center_A == tnsr::I<double, 3, Frame::Grid>{{8.0, 0.0, 0.0}});
  CHECK(grid_center_B == tnsr::I<double, 3, Frame::Grid>{{-5.5, 0.0, 0.0}});

  const std::unique_ptr<DomainCreator<3>> creator_no_excision =
      std::make_unique<domain::creators::Brick>(
          std::array{0.0, 0.0, 0.0}, std::array{1.0, 1.0, 1.0},
          std::array{0_st, 0_st, 0_st}, std::array{2_st, 2_st, 2_st},
          std::array{false, false, false});

  CHECK_THROWS_WITH(
      Tags::ExcisionCenter<ObjectLabel::B>::create_from_options(
          creator_no_excision),
      Catch::Contains(" is not in the domains excision spheres but is needed "
                      "to generate the ExcisionCenter"));
}

template <bool Override>
struct Metavariables {
  static constexpr size_t volume_dim = 3;
  static constexpr bool override_functions_of_time = Override;
};

void test_functions_of_time() {
  TestHelpers::db::test_simple_tag<Tags::FunctionsOfTimeInitialize>(
      "FunctionsOfTime");

  CHECK(std::is_same_v<
        Tags::FunctionsOfTimeInitialize::option_tags<Metavariables<true>>,
        tmpl::list<
            domain::OptionTags::DomainCreator<Metavariables<true>::volume_dim>,
            domain::FunctionsOfTime::OptionTags::FunctionOfTimeFile,
            domain::FunctionsOfTime::OptionTags::FunctionOfTimeNameMap>>);

  CHECK(std::is_same_v<
        Tags::FunctionsOfTimeInitialize::option_tags<Metavariables<false>>,
        tmpl::list<domain::OptionTags::DomainCreator<
            Metavariables<false>::volume_dim>>>);
}

void test_size_shape() {
  domain::FunctionsOfTime::register_derived_with_charm();

  using ShapeTrans = domain::CoordinateMaps::ShapeMapTransitionFunctions::
      ShapeMapTransitionFunction;
  using SphereTrans =
      domain::CoordinateMaps::ShapeMapTransitionFunctions::SphereTransition;
  using Shape = domain::CoordinateMaps::TimeDependent::Shape;
  using Size =
      domain::CoordinateMaps::TimeDependent::SphericalCompression<false>;
  using FoT = domain::FunctionsOfTime::FunctionOfTime;
  using PP = domain::FunctionsOfTime::PiecewisePolynomial<1>;

  const double initial_time = 0.0;
  const double expr_time = 10.0;

  const double inner_radius = 1.5;
  const double outer_radius = 2.0;

  const std::array<double, 3> center{0.0, 0.0, 0.0};
  const size_t l_max = 10;
  std::unique_ptr<ShapeTrans> shape_trans =
      std::make_unique<SphereTrans>(inner_radius, outer_radius);

  Shape shape_map{center, l_max, l_max, std::move(shape_trans), "Shape"};
  Size size_map{"Size", inner_radius, outer_radius, center};

  DataVector shape_coefs{YlmSpherepack::spectral_size(l_max, l_max), 0.0};
  DataVector size_coefs{1, 0.0};
  SpherepackIterator iter{l_max, l_max};

  const double initial_size = 4.0;

  shape_coefs[iter.set(0, 0)()] =
      sqrt(2.0 / M_PI) * initial_size / inner_radius;
  size_coefs[0] = initial_size;

  std::unordered_map<std::string, std::unique_ptr<FoT>> functions_of_time{};

  functions_of_time["Shape"] = std::make_unique<PP>(
      initial_time,
      std::array<DataVector, 2>{
          {shape_coefs, DataVector{shape_coefs.size(), 0.0}}},
      expr_time);
  functions_of_time["Size"] = std::make_unique<PP>(
      initial_time, std::array<DataVector, 2>{{size_coefs, DataVector{1, 0.0}}},
      expr_time);

  const std::array<double, 3> input_coord{1.6, 0.0, 0.0};

  const std::array<double, 3> size_output_coord =
      size_map(input_coord, initial_time, functions_of_time);
  const std::array<double, 3> shape_output_coord =
      shape_map(input_coord, initial_time, functions_of_time);

  CHECK(shape_output_coord == size_output_coord);
}

SPECTRE_TEST_CASE("Unit.Domain.Creators.Tags", "[Unit][Domain]") {
  test_simple_tags<1>();
  test_simple_tags<2>();
  test_simple_tags<3>();

  test_center_tags();

  test_functions_of_time();

  test_size_shape();
}
}  // namespace
}  // namespace domain
