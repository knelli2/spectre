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
#include "Framework/ActionTesting.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/IO/Observers/MockH5.hpp"
#include "Helpers/IO/Observers/MockWriteReductionDataRow.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/RegisterDerivedWithCharm.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/Sphere.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
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
using tags_to_observe_on_target = tmpl::list<
    gr::surfaces::Tags::SurfaceIntegralCompute<Tags::Negate, Frame::NoFrame>,
    gr::surfaces::Tags::SurfaceIntegralCompute<Tags::Square, Frame::NoFrame>>;
using non_observation_tags =
    tmpl::list<Tags::NegateCompute, Tags::SquareCompute>;
using volume_compute_tags = tmpl::list<Tags::TestSolution>;

struct Target;

using Callback = intrp2::callbacks::ObserveTimeSeriesOnSurface<
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

struct TestMetavars {
  using component_list =
      tmpl::list<TestHelpers::observers::MockObserverWriter<TestMetavars>>;

  struct factory_creation {
    using factory_classes = tmpl::map<
        tmpl::pair<intrp2::callbacks::Callback<Target>, tmpl::list<Callback>>>;
  };
};

void test_option_parsing() {
  intrp2::callbacks::register_derived_with_charm<tmpl::list<Target>>();
  const auto callback = TestHelpers::test_creation<Callback>(
      "SubfileName: Crab\n"
      "ValuesToObserve:\n"
      "  - SurfaceIntegral(Negate)\n"
      "  - SurfaceIntegral(Square)");

  const std::unordered_set<std::string> expected_observables{
      "SurfaceIntegral(Square)", "SurfaceIntegral(Negate)"};
  CHECK(expected_observables == callback.observables());

  CHECK_THROWS_WITH(
      TestHelpers::test_creation<Callback>("SubfileName: Crab\n"
                                           "ValuesToObserve:\n"
                                           "  - SurfaceIntegral(Negate)\n"
                                           "  - Square"),
      Catch::Matchers::ContainsSubstring("Invalid selection:") and
          Catch::Matchers::ContainsSubstring("Possible choices are:"));
  CHECK_THROWS_WITH(
      TestHelpers::test_creation<Callback>("SubfileName: Crab\n"
                                           "ValuesToObserve:\n"
                                           "  - SurfaceIntegral(Negate)\n"
                                           "  - SurfaceIntegral(Negate)"),
      Catch::Matchers::ContainsSubstring("specified multiple times"));

  const auto callback_ptr =
      TestHelpers::test_factory_creation<intrp2::callbacks::Callback<Target>,
                                         Callback>(
          "ObserveTimeSeriesOnSurface:\n"
          "  SubfileName: Crab\n"
          "  ValuesToObserve:\n"
          "    - SurfaceIntegral(Negate)\n"
          "    - SurfaceIntegral(Square)");
  CHECK(expected_observables == callback_ptr->observables());
  const auto serialized_callback_ptr = serialize_and_deserialize(callback_ptr);
  CHECK(expected_observables == serialized_callback_ptr->observables());

  const auto cloned_ptr = callback_ptr->get_clone();
  CHECK(expected_observables == cloned_ptr->observables());
}

void test_callback() {
  using observer_writer =
      TestHelpers::observers::MockObserverWriter<TestMetavars>;
  const auto callback_ptr =
      TestHelpers::test_factory_creation<intrp2::callbacks::Callback<Target>,
                                         Callback>(
          "ObserveTimeSeriesOnSurface:\n"
          "  SubfileName: Crab\n"
          "  ValuesToObserve:\n"
          "    - SurfaceIntegral(Negate)\n"
          "    - SurfaceIntegral(Square)");

  ActionTesting::MockRuntimeSystem<TestMetavars> runner{{}};
  ActionTesting::set_phase(make_not_null(&runner),
                           Parallel::Phase::Initialization);
  ActionTesting::emplace_nodegroup_component_and_initialize<observer_writer>(
      make_not_null(&runner), {});

  auto& cache = ActionTesting::cache<observer_writer>(runner, 0_st);

  // We aren't testing the compute tag structure here, so we only add the simple
  // tags and just put random values.
  auto intrp_box = db::create<tmpl::list<
      gr::surfaces::Tags::SurfaceIntegral<Tags::Negate, ::Frame::NoFrame>,
      gr::surfaces::Tags::SurfaceIntegral<Tags::Square, ::Frame::NoFrame>>>(
      -12345.0, 54321.0);
  const db::Access& intrp_access = db::as_access(intrp_box);

  const double time = 1.5;
  callback_ptr->invoke(intrp_access, cache, time);

  CHECK(ActionTesting::number_of_queued_threaded_actions<observer_writer>(
            runner, 0_st) == 1);
  ActionTesting::invoke_queued_threaded_action<observer_writer>(
      make_not_null(&runner), 0_st);

  const TestHelpers::observers::MockH5File& mock_h5 =
      ActionTesting::get_databox_tag<
          observer_writer, TestHelpers::observers::MockReductionFileTag>(runner,
                                                                         0_st);
  const TestHelpers::observers::MockDat& mock_dat = mock_h5.get_dat("/Crab");
  const Matrix& mock_data = mock_dat.get_data();

  CHECK(mock_dat.get_legend() ==
        std::vector<std::string>{"Time", "SurfaceIntegral(Negate)",
                                 "SurfaceIntegral(Square)"});
  CHECK(mock_data.columns() == 3);
  CHECK(mock_data.rows() == 1);
  CHECK(mock_data(0, 0) == time);
  CHECK(mock_data(0, 1) == -12345.0);
  CHECK(mock_data(0, 2) == 54321.0);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.Callbacks."
    "ObserveTimeSeriesOnSurface",
    "[Unit]") {
  test_option_parsing();
  test_callback();
}
