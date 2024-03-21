// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>
#include <pup.h>
#include <string>
#include <unordered_set>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Framework/ActionTesting.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Actions/InterpolationTargetVarsFromElement.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/LineSegment.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace {
namespace Tags {
struct Time : db::SimpleTag {
  using type = double;
};
}  // namespace Tags

template <typename Target>
struct MockCallback : public intrp2::callbacks::Callback<Target>,
                      public tt::ConformsTo<intrp2::protocols::Callback> {
  using tags_to_observe_on_target = tmpl::list<>;
  using non_observation_tags_on_target = tmpl::list<>;
  using volume_compute_tags = tmpl::list<>;

  MockCallback() = default;
  explicit MockCallback(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(intrp2::callbacks::Callback<Target>,
                                     MockCallback);

  using options = tmpl::list<>;
  static constexpr Options::String help = {"halp"};

  void pup(PUP::er& p) override { p | values_to_observe_; }

  const std::unordered_set<std::string>& observables() const override {
    return values_to_observe_;
  }

  template <typename Metavariables>
  void apply(const db::Access& access,
             Parallel::GlobalCache<Metavariables>& cache, const double time) {
    (void)access;
    (void)cache;
    (void)time;
  }

 private:
  std::unordered_set<std::string> values_to_observe_{};
};

template <typename Target>
// NOLINTNEXTLINE
PUP::able::PUP_ID MockCallback<Target>::my_PUP_ID = 0;

using Points = intrp2::points::LineSegment<1>;

struct MockTarget : tt::ConformsTo<intrp2::protocols::Target> {
  using temporal_id_tag = Tags::Time;
  using frame = NoSuchType;
  using points = Points;
  using possible_runtime_callbacks = tmpl::list<MockCallback<MockTarget>>;
  using compile_time_callbacks = tmpl::list<>;
};

template <typename Metavariables>
struct MockTargetComponent {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = std::string;
  using const_global_cache_tags = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 //  tmpl::list<intrp2::Tags::DbAccess>>>>>;
                                 tmpl::list<>>>>>;
};

struct TestMetavars {
  static constexpr size_t volume_dim = 1;

  using component_list = tmpl::list<MockTargetComponent<TestMetavars>>;
};

void test() {
  using target_component = MockTargetComponent<TestMetavars>;
  (void)TestMetavars::volume_dim;

  using tag_list = tmpl::list<>;
  using BoxType = db::compute_databox_type<tag_list>;
  std::unique_ptr<db::Access> access = std::make_unique<BoxType>();

  ActionTesting::MockRuntimeSystem<TestMetavars> runner{{}};
  // ActionTesting::emplace_array_component_and_initialize<target_component>(
  //     make_not_null(&runner), ActionTesting::NodeId{0},
  //     ActionTesting::LocalCoreId{0}, "MockTarget", {std::move(access)});
  ActionTesting::emplace_array_component<target_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, "MockTarget");
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.VarsFromElement", "[Unit]") {
  test();
}
