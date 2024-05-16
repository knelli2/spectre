// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>
#include <unordered_set>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/AsAccess.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/IO/Observers/MockH5.hpp"
#include "Helpers/IO/Observers/MockWriteReductionDataRow.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/TensorData.hpp"
#include "IO/Observer/Initialize.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/StrahlkorperFunctions.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/ObserveSurfaceData.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/RegisterDerivedWithCharm.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/Sphere.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
namespace Tags {
struct Time : db::SimpleTag {
  using type = double;
};
struct TestSolution : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct Square : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct SquareCompute : Square, db::ComputeTag {
  static void function(gsl::not_null<Scalar<DataVector>*> result,
                       const Scalar<DataVector>& x) {
    get(*result) = square(get(x));
  }
  using argument_tags = tmpl::list<TestSolution>;
  using base = Square;
  using return_type = Scalar<DataVector>;
};
struct Negate : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct NegateCompute : Negate, db::ComputeTag {
  static void function(gsl::not_null<Scalar<DataVector>*> result,
                       const Scalar<DataVector>& x) {
    get(*result) = -get(x);
  }
  using argument_tags = tmpl::list<Square>;
  using base = Negate;
  using return_type = Scalar<DataVector>;
};
}  // namespace Tags

using tags_on_target_from_points = tmpl::list<>;
using tags_to_observe_on_target = tmpl::list<Tags::Negate, Tags::Square>;
using non_observation_tags =
    tmpl::list<Tags::NegateCompute, Tags::SquareCompute>;
using volume_compute_tags = tmpl::list<Tags::TestSolution>;

struct Target;

using Callback = intrp2::callbacks::ObserveDataOnStrahlkorper<
    Target, tags_to_observe_on_target, non_observation_tags,
    volume_compute_tags>;
static_assert(tt::assert_conforms_to_v<Callback, intrp2::protocols::Callback>);

struct Target {
  using temporal_id_tag = Tags::Time;
  using frame = NoSuchType;
  using points = intrp2::points::Sphere;
  using possible_runtime_callbacks = tmpl::list<Callback>;
  using compile_time_callbacks = tmpl::list<>;
};

struct MockWriteVolumeData {
  template <typename ParallelComponent, typename DbTagsList,
            typename Metavariables, typename ArrayIndex, typename... Ts,
            typename DataBox = db::DataBox<DbTagsList>>
  static void apply(db::DataBox<DbTagsList>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const gsl::not_null<Parallel::NodeLock*> /*node_lock*/,
                    const std::string& /*h5_file_name*/,
                    const std::string& subfile_path,
                    const observers::ObservationId& /*observation_id*/,
                    std::vector<ElementVolumeData>&& volume_data) {
    CHECK(subfile_path == "/Lobster");
    CHECK(volume_data.size() == 1);
    CHECK(volume_data[0].element_name == "Lobster");
    CHECK(volume_data[0].basis ==
          std::vector<Spectral::Basis>{2, Spectral::Basis::SphericalHarmonic});
    CHECK(volume_data[0].quadrature ==
          std::vector<Spectral::Quadrature>{Spectral::Quadrature::Gauss,
                                            Spectral::Quadrature::Equiangular});
    // 2 for the scalars, 3 for inertial coords, then possible another 3 for the
    // other frame coords
    CHECK((volume_data[0].tensor_components.size() == 5 or
           volume_data[0].tensor_components.size() == 8));
    std::unordered_set<std::string> expected_names{"Negate", "Square"};
    for (const auto& tensor_component : volume_data[0].tensor_components) {
      CHECK((expected_names.contains(tensor_component.name) or
             tensor_component.name.find("Coordinates") != std::string::npos));
    }
  }
};

template <typename Metavariables>
struct MockObserverWriter {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockNodeGroupChare;
  using array_index = size_t;
  using const_global_cache_tags = tmpl::list<observers::Tags::ReductionFileName,
                                             observers::Tags::SurfaceFileName>;
  using simple_tags =
      typename observers::Actions::InitializeWriter<Metavariables>::simple_tags;
  using compute_tags = typename observers::Actions::InitializeWriter<
      Metavariables>::compute_tags;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<tmpl::list<
                         TestHelpers::observers::MockReductionFileTag>>,
                     observers::Actions::InitializeWriter<Metavariables>>>,
      Parallel::PhaseActions<Parallel::Phase::Testing, tmpl::list<>>>;

  using component_being_mocked = observers::ObserverWriter<Metavariables>;

  using replace_these_threaded_actions =
      tmpl::list<::observers::ThreadedActions::WriteReductionDataRow,
                 ::observers::ThreadedActions::WriteVolumeData>;
  using with_these_threaded_actions =
      tmpl::list<TestHelpers::observers::MockWriteReductionDataRow,
                 MockWriteVolumeData>;
};

struct TestMetavars {
  using component_list = tmpl::list<MockObserverWriter<TestMetavars>>;

  using const_global_cache_tags = tmpl::list<::domain::Tags::Domain<3>>;
  using mutable_global_cache_tags =
      tmpl::list<::domain::Tags::FunctionsOfTimeInitialize>;

  using observed_reduction_data_tags = tmpl::list<>;

  struct factory_creation {
    using factory_classes = tmpl::map<
        tmpl::pair<intrp2::callbacks::Callback<Target>, tmpl::list<Callback>>>;
  };
};

using observer_writer = MockObserverWriter<TestMetavars>;

void test_option_parsing() {
  intrp2::callbacks::register_derived_with_charm<tmpl::list<Target>>();
  const auto callback = TestHelpers::test_creation<Callback>(
      "SubfileName: Lobster\n"
      "ValuesToObserve:\n"
      "  - Negate\n"
      "  - Square");

  const std::unordered_set<std::string> expected_observables{"Square",
                                                             "Negate"};
  CHECK(expected_observables == callback.observables());

  CHECK_THROWS_WITH(
      TestHelpers::test_creation<Callback>("SubfileName: Lobster\n"
                                           "ValuesToObserve:\n"
                                           "  - SurfaceIntegral(Negate)\n"
                                           "  - Square"),
      Catch::Matchers::ContainsSubstring("Invalid selection:") and
          Catch::Matchers::ContainsSubstring("Possible choices are:"));
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<Callback>("SubfileName: Lobster\n"
                                           "ValuesToObserve:\n"
                                           "  - Negate\n"
                                           "  - Negate"),
      Catch::Matchers::ContainsSubstring("specified multiple times"));

  const auto callback_ptr =
      TestHelpers::test_factory_creation<intrp2::callbacks::Callback<Target>,
                                         Callback>(
          "ObserveDataOnStrahlkorper:\n"
          "  SubfileName: Lobster\n"
          "  ValuesToObserve:\n"
          "    - Negate\n"
          "    - Square");
  CHECK(expected_observables == callback_ptr->observables());
  const auto serialized_callback_ptr = serialize_and_deserialize(callback_ptr);
  CHECK(expected_observables == serialized_callback_ptr->observables());

  const auto cloned_ptr = callback_ptr->get_clone();
  CHECK(expected_observables == cloned_ptr->observables());
}

template <typename Runner>
void check_ylm_coefs(const gsl::not_null<Runner*> runner,
                     const ylm::Strahlkorper<Frame::NoFrame>& expected_surface,
                     const std::string& frame, const double time) {
  // Parameters chosen to match SurfaceE choices below
  const size_t l_max = expected_surface.l_max();
  const std::array<double, 3>& expansion_center =
      expected_surface.expansion_center();

  const std::vector<std::string> ylm_expected_legend{
      "Time",
      frame + "ExpansionCenter_x",
      frame + "ExpansionCenter_y",
      frame + "ExpansionCenter_z",
      "Lmax",
      "coef(0,0)",
      "coef(1,-1)",
      "coef(1,0)",
      "coef(1,1)",
      "coef(2,-2)",
      "coef(2,-1)",
      "coef(2,0)",
      "coef(2,1)",
      "coef(2,2)",
      "coef(3,-3)",
      "coef(3,-2)",
      "coef(3,-1)",
      "coef(3,0)",
      "coef(3,1)",
      "coef(3,2)",
      "coef(3,3)"};
  const size_t expected_num_columns = ylm_expected_legend.size();

  // Check that the H5 file was written correctly.
  const auto& file = ActionTesting::get_databox_tag<
      observer_writer, TestHelpers::observers::MockReductionFileTag>(*runner,
                                                                     0);
  const auto& ylm_dat_file = file.get_dat("/Lobster_Ylm");
  const Matrix& ylm_written_data = ylm_dat_file.get_data();
  const auto& ylm_written_legend = ylm_dat_file.get_legend();

  CHECK(ylm_written_legend.size() == expected_num_columns);
  CHECK(ylm_written_data.columns() == expected_num_columns);
  CHECK(ylm_written_legend == ylm_expected_legend);

  std::vector<double> ylm_expected_data{
      time, expansion_center[0], expansion_center[1], expansion_center[2],
      static_cast<double>(l_max)};

  ylm::SpherepackIterator iter(l_max, l_max);
  for (size_t l = 0; l <= l_max; l++) {
    for (int m = -static_cast<int>(l); m <= static_cast<int>(l); m++) {
      iter.set(l, m);
      ylm_expected_data.push_back(expected_surface.coefficients()[iter()]);
    }
  }

  for (size_t i = 0; i < expected_num_columns; i++) {
    CHECK(ylm_written_data(0, i) == ylm_expected_data[i]);
  }
}

void test_callback(const std::string& frame) {
  domain::creators::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();
  const std::string reduction_filename{"ReductionJonBonJovi"};
  const std::string surface_filename{"SurfaceJonBonJovi"};

  const double time = 1.5;
  // No actual time dependent maps here. This is ok because the function
  // strahlkorper_coords_in_different_frame can handle time independent blocks
  // and we're testing that the callback is called correctly, not that the frame
  // transformations work. Those are tested elsewhere.
  const domain::creators::Sphere sphere{
      0.1, 10.0, domain::creators::Sphere::Excision{}, 0_st, 5_st, false};

  const auto callback_ptr =
      TestHelpers::test_factory_creation<intrp2::callbacks::Callback<Target>,
                                         Callback>(
          "ObserveDataOnStrahlkorper:\n"
          "  SubfileName: Lobster\n"
          "  ValuesToObserve:\n"
          "    - Negate\n"
          "    - Square");

  ActionTesting::MockRuntimeSystem<TestMetavars> runner{
      {sphere.create_domain(), reduction_filename, surface_filename},
      {sphere.functions_of_time()}};
  ActionTesting::set_phase(make_not_null(&runner),
                           Parallel::Phase::Initialization);
  ActionTesting::emplace_nodegroup_component_and_initialize<observer_writer>(
      make_not_null(&runner), {});

  auto& cache = ActionTesting::cache<observer_writer>(runner, 0_st);

  const size_t l_max = 3;
  const double radius = 1.2;
  const std::array<double, 3> center{1.0, 2.0, 3.0};
  ylm::Strahlkorper<Frame::NoFrame> strahlkorper{l_max, radius, center};
  tnsr::I<DataVector, 3, Frame::NoFrame> points =
      ylm::cartesian_coords(strahlkorper, ylm::radius(strahlkorper),
                            ylm::rhat(ylm::theta_phi(strahlkorper)));

  // We aren't testing the compute tag structure here, so we only add the
  // simple tags and just put some values.
  auto intrp_box =
      db::create<tmpl::list<intrp2::Tags::Frame, intrp2::Tags::Points<3>,
                            ylm::Tags::Strahlkorper<Frame::NoFrame>,
                            Tags::Negate, Tags::Square>>(
          frame, points, strahlkorper,
          Scalar<DataVector>{get<0>(points).size()},
          Scalar<DataVector>{get<0>(points).size()});
  const db::Access& intrp_access = db::as_access(intrp_box);

  (void)cache;
  (void)intrp_access;

  callback_ptr->invoke(intrp_access, cache, time);

  const size_t number_queued_actions =
      ActionTesting::number_of_queued_threaded_actions<observer_writer>(runner,
                                                                        0);
  CHECK(number_queued_actions == 2);

  for (size_t i = 0; i < number_queued_actions; i++) {
    ActionTesting::invoke_queued_threaded_action<observer_writer>(
        make_not_null(&runner), 0);
  }

  check_ylm_coefs(make_not_null(&runner), strahlkorper, frame, time);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.Callbacks."
    "ObserveDataOnSurface",
    "[Unit]") {
  test_option_parsing();
  for (const auto& frame : {"Grid", "Distorted", "Inertial"}) {
    test_callback(frame);
  }
}
