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

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Domain/Structure/CreateInitialMesh.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/Runtime/PointsHelpers.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/Sphere.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/CartesianProduct.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Spherepack.hpp"
#include "Utilities/TMPL.hpp"

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
        const double phi = two_pi_over_n_phi * static_cast<double>(i_phi);
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
          const double phi = two_pi_over_n_phi * static_cast<double>(i_phi);
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
  CHECK(sphere.number_of_sets_of_points() == number_of_spheres);
}

struct TestMetavars {
  using component_list = tmpl::list<>;

  using const_global_cache_tags = tmpl::list<domain::Tags::Domain<3>>;
};

void test_operator() {
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();

  const double time = 0.0;
  // We do this so we know exactly where some points are
  const size_t l_max = 2;
  // There are (l_max + 1)*(2l_max + 1) = 15 collocation points
  const size_t num_points = ylm::Spherepack::physical_size(l_max, l_max);
  (void)num_points;
  const std::array<double, 3> center{0.0, 0.0, 0.0};
  const std::set<double> radii_set{1.0, 2.0};
  const intrp2::AngularOrdering angular_ordering =
      intrp2::AngularOrdering::Strahlkorper;
  const intrp2::points::Sphere sphere{
      l_max, center, std::vector<double>{radii_set.begin(), radii_set.end()},
      angular_ordering};

  // One radius is inside the domain, one is out
  domain::creators::Sphere sphere_creator{
      1.5, 10.0, domain::creators::Sphere::Excision{}, 0_st, 4_st, true};
  Parallel::GlobalCache<TestMetavars> cache{{sphere_creator.create_domain()}};
  const auto domain = sphere_creator.create_domain();
  const auto initial_extents = sphere_creator.initial_extents();
  std::vector<ElementId<3>> element_ids =
      initial_element_ids(sphere_creator.initial_refinement_levels());

  using BoxType = db::compute_databox_type<tmpl::list<
      domain::Tags::Mesh<3>, domain::Tags::ElementMap<3, Frame::Grid>,
      domain::Tags::ElementMap<3, Frame::Distorted>,
      domain::Tags::ElementMap<3, Frame::Inertial>>>;
  BoxType box{};

  // We don't have time dependent maps so we just use the Grid frame
  const auto expected_block_logical_coordinates =
      block_logical_coordinates_in_frame(
          cache, time, sphere.target_points_no_frame(), "Grid");
  // Blocks in order are +/-z, +/-y, +/-x
  // We only have three points in the theta direction, so +/-z blocks have all
  // the endpoints (5 since we have 5 phi points). Then just based of the phi
  // spacing of the 5 middle-theta points, +/-y and +x have 1 point, and -x has
  // two points.
  const std::vector<size_t> expected_block_coord_sizes{5, 5, 1, 1, 1, 2};

  // Same number of elements as blocks
  for (size_t block = 0; block < element_ids.size(); block++) {
    CAPTURE(block);
    const ElementId<3>& element_id = element_ids[block];
    const Element<3> element{element_id, {}};
    db::mutate<domain::Tags::Mesh<3>, domain::Tags::ElementMap<3, Frame::Grid>,
               domain::Tags::ElementMap<3, Frame::Distorted>,
               domain::Tags::ElementMap<3, Frame::Inertial>>(
        [&](const gsl::not_null<Mesh<3>*> mesh,
            const gsl::not_null<ElementMap<3, Frame::Grid>*> element_map_grid,
            const gsl::not_null<ElementMap<3, Frame::Distorted>*>
                element_map_distorted,
            const gsl::not_null<ElementMap<3, Frame::Inertial>*>
                element_map_inertial) {
          // This is a comment
          *mesh = domain::Initialization::create_initial_mesh(
              initial_extents, element.id(),
              Spectral::Quadrature::GaussLobatto);
          *element_map_grid =
              ElementMap<3, Frame::Grid>{element.id(), domain.blocks()[block]};
          *element_map_distorted = ElementMap<3, Frame::Distorted>{
              element.id(), domain.blocks()[block]};
          *element_map_inertial = ElementMap<3, Frame::Inertial>{
              element.id(), domain.blocks()[block]};
        },
        make_not_null(&box));
    const auto test_frame = [&]<typename LocalFrame>() {
      const std::string frame = pretty_type::name<LocalFrame>();
      INFO(frame + " frame");
      ObservationBox<
          tmpl::list<domain::Tags::LogicalCoordinates<3>,
                     domain::Tags::MappedCoordinates<
                         domain::Tags::ElementMap<3, LocalFrame>,
                         domain::Tags::Coordinates<3, Frame::ElementLogical>>>,
          BoxType>
          obs_box{make_not_null(&box)};

      const auto block_coords = sphere(cache, time, frame, element, obs_box);
      CHECK(block_coords.has_value());
      CHECK(block_coords->size() == expected_block_coord_sizes[block]);
      const domain::BlockId block_id{block};
      for (const auto& block_coord : block_coords.value()) {
        CHECK(block_coord.has_value());
        CHECK(block_coord->id == block_id);
        CHECK(alg::any_of(
            expected_block_logical_coordinates,
            [&block_coord, &block_id](const auto& expected_block_coord) {
              // Make sure we are checking the correct block first before we
              // check the logical coordinate
              return expected_block_coord.has_value() and
                     expected_block_coord->id == block_id and
                     expected_block_coord->data == block_coord->data;
            }));
      }

      {
        INFO(frame + " frame error");
        ObservationBox<tmpl::list<>, BoxType> bad_obs_box{make_not_null(&box)};

        CHECK_THROWS_WITH(sphere(cache, time, frame, element, bad_obs_box),
                          Catch::Matchers::ContainsSubstring(
                              "Could not find a coordinates tag "
                              "in the box for the " +
                              frame));
      }
    };

    // We aren't testing that block_logical_coords works, just that the sphere
    // has the correct behavior with the different frames. So for this test all
    // points are the same for each frame
    test_frame.operator()<Frame::Grid>();
    test_frame.operator()<Frame::Distorted>();
    test_frame.operator()<Frame::Inertial>();
  }
}

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
  test_operator();
}
