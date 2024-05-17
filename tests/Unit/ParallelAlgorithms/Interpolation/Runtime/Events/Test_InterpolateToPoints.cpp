
// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/AccessWrapper.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Actions/InterpolationTargetVarsFromElement.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Events/InterpolateToPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/Sphere.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Surfaces/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/NoSuchType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

#include "Parallel/Printf/Printf.hpp"

namespace {
namespace Tags {
struct Time : db::SimpleTag {
  using type = double;
};
}  // namespace Tags

struct Target;

// specialize to double for now
template <typename Tag>
struct Dot : db::SimpleTag {
  static std::string name() { return "Dot("s + db::tag_name<Tag>() + ")"s; }
  using type = double;
};

template <typename Tag>
struct DotCompute : Dot<Tag>, db::ComputeTag {
  using base = Dot<Tag>;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<Tag>;

  static void function(const gsl::not_null<return_type*> result,
                       const typename Tag::type& tag) {
    *result = dot(tag);
  }
};

using tags_to_observe =
    tmpl::list<DotCompute<gr::Tags::Shift<double, 3>>,
               ylm::Tags::EuclideanSurfaceIntegralCompute<
                   Dot<gr::Tags::Lapse<double>>, Frame::NoFrame>>;
using non_observation_tags = tmpl::list<>;
using volume_compute_tags = tmpl::list<>;
using TestCallback = intrp2::callbacks::ObserveTimeSeriesOnSurface<
    Target, tags_to_observe, non_observation_tags, volume_compute_tags>;

struct Target : public tt::ConformsTo<intrp2::protocols::Target> {
  using temporal_id_tag = Tags::Time;

  using frame = NoSuchType;
  using points = intrp2::points::Sphere;

  using possible_runtime_callbacks = tmpl::list<TestCallback>;
  using compile_time_callbacks = tmpl::list<TestCallback>;
};

using Event = intrp2::Events::InterpolateToPoints<Target, 3>;

// We are only testing the Event
struct MockReceiveVolumeTensors {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename TemporalId>
  static void apply(
      db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& cache,
      const std::string& array_index, const TemporalId& temporal_id,
      const std::unordered_map<std::string, std::vector<DataVector>>&
          received_interpolated_vars,
      const std::vector<size_t>& global_offsets,
      const bool check_for_invalid_points,
      const bool vars_have_already_been_received = false) {
    (void)box;
    (void)cache;
    (void)array_index;
    (void)temporal_id;
    (void)received_interpolated_vars;
    (void)global_offsets;
    (void)check_for_invalid_points;
    (void)vars_have_already_been_received;
  }
};

template <typename Metavariables>
struct MockTargetComponent {
  using metavariables = Metavariables;
  using component_being_mocked = intrp2::InterpolationTargets<metavariables>;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = std::string;
  using const_global_cache_tags = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<tmpl::list<>>>>>;

  using replace_these_simple_actions =
      tmpl::list<intrp2::Actions::ReceiveVolumeTensors<Target>>;
  using with_these_simple_actions = tmpl::list<MockReceiveVolumeTensors>;
};

struct TestMetavars {
  static constexpr size_t volume_dim = 3;

  using component_list = tmpl::list<MockTargetComponent<TestMetavars>>;
};

void test() {
  constexpr size_t Dim = TestMetavars::volume_dim;
  ActionTesting::MockRuntimeSystem<TestMetavars> runner{{}};

  const size_t l_max = 4;
  const std::array center{0.0, 0.0, 0.0};
  const std::vector<double> radii{1.0, 2.0};
  const intrp2::AngularOrdering angular_ordering =
      intrp2::AngularOrdering::Strahlkorper;
  intrp2::points::Sphere sphere{l_max, center, radii, angular_ordering};

  const double invalid_fill_value = -1.0;
  const std::string frame{"Inertial"};

  const std::string subfile{"subfile"};
  const std::vector<std::string> vars_to_observe{"Dot(Shift)"};
  std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>>
      expected_callbacks(2);
  std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>> callbacks(
      2);

  callbacks[0] = std::make_unique<TestCallback>(subfile, vars_to_observe);
  callbacks[1] = std::make_unique<TestCallback>(subfile, vars_to_observe);
  expected_callbacks[0] =
      std::make_unique<TestCallback>(subfile, vars_to_observe);
  expected_callbacks[1] =
      std::make_unique<TestCallback>(subfile, vars_to_observe);

  intrp2::Events::InterpolateToPoints<Target, Dim> event{
      sphere, invalid_fill_value, frame, std::move(callbacks)};

  {
    INFO("Check callback observables");
    std::unique_ptr<intrp2::AccessWrapper> access_wrapper =
        std::make_unique<intrp2::TargetAccessWrapper<Target, Dim>>();
    event.initialize_target_element_box(make_not_null(&access_wrapper));

    db::Access& access = access_wrapper->access();

    CHECK(db::get<intrp2::Tags::Frame>(access) == frame);
    CHECK(db::get<intrp2::Tags::Points<Dim>>(access) ==
          sphere.target_points_no_frame());
    CHECK(db::get<intrp2::Tags::InvalidPointsFillValue>(access) ==
          invalid_fill_value);
    const auto& access_callbacks =
        db::get<intrp2::Tags::Callbacks<Target>>(access);
    CHECK(access_callbacks.size() == expected_callbacks.size());
    for (size_t i = 0; i < expected_callbacks.size(); i++) {
      CHECK(access_callbacks[i]->observables() ==
            expected_callbacks[i]->observables());
    }
  }
  {
    INFO("Check tag mapping");
    const auto& tensors_to_observe = event.tensors_to_observe();
    const auto& volume_tensors_to_send = event.volume_tensors_to_send();
    const auto& all_points_tags = event.all_points_tags();
    const auto& non_obs_tags_to_volume_tags =
        event.non_obs_tags_to_volume_tags();
    const auto& observation_tags_to_volume_tags =
        event.observation_tags_to_volume_tags();

    Parallel::printf(
        "tensors_to_observe: %s\n"
        "volume_tensors_to_send: %s\n"
        "all_points_tags: %s\n"
        "non_obs_tags_to_volume_tags: %s\n"
        "observation_tags_to_volume_tags: %s\n",
        tensors_to_observe, volume_tensors_to_send, all_points_tags,
        non_obs_tags_to_volume_tags, observation_tags_to_volume_tags);

    CHECK(alg::all_of(vars_to_observe, [&](const auto& var) {
      return tensors_to_observe.contains(var);
    }));
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.InterpolateToPoints",
    "[Unit]") {
  test();
}
