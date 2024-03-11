// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <deque>
#include <memory>
#include <pup.h>
#include <string>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Variables.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Events/MarkAsInterpolation.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Metavariables.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

class DataVector;

namespace {
struct AccessTag : db::SimpleTag {
  using type = bool;
};

template <size_t Index>
class TestEvent : public ::Event, public intrp2::Events::MarkAsInterpolation {
 public:
  explicit TestEvent(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(TestEvent);  // NOLINT

  std::string name() const override {
    return "TestEvent" + std::to_string(Index);
  }

  TestEvent() = default;

  std::unique_ptr<db::Access> initialize_target_element_box() const override {
    auto box = db::create<tmpl::list<AccessTag>>(true);
    return std::make_unique<db::DataBox<tmpl::list<AccessTag>>>(std::move(box));
  }

  bool needs_evolved_variables() const override { return false; }

  // NOLINTNEXTLINE
  void pup(PUP::er& /*p*/) override {}
};

template <size_t Index>
PUP::able::PUP_ID TestEvent<Index>::my_PUP_ID = 0;  // NOLINT

template <typename Metavariables>
struct mock_interpolation_target {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = std::string;
  using const_global_cache_tags = tmpl::list<::Tags::EventsAndTriggers>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<intrp2::Actions::InitializeInterpolationTarget>>>;
};

constexpr size_t num_elements = 3;
template <typename SizeType>
struct get_event {
  using type = TestEvent<SizeType::value>;
};
using event_list =
    tmpl::transform<tmpl::range<size_t, 0, num_elements>, get_event<tmpl::_1>>;

struct Allocator : tt::ConformsTo<Parallel::protocols::ArrayElementsAllocator> {
  template <typename Component>
  using array_allocation_tags = tmpl::list<>;
  template <typename ParallelComponent, typename Metavariables,
            typename... InitializationTags>
  static void apply(
      Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/,
      const tuples::TaggedTuple<
          InitializationTags...>& /*initialization_items*/,
      const tuples::tagged_tuple_from_typelist<
          typename ParallelComponent::array_allocation_tags>&
      /*array_allocation_items*/
      = {},
      const std::unordered_set<size_t>& /*procs_to_ignore*/ = {}) {}
};

struct Metavariables {
  using component_list = tmpl::list<mock_interpolation_target<Metavariables>>;

  struct intrp : tt::ConformsTo<intrp2::protocols::Metavariables> {
    // We don't need the allocator for this test so we just use one that doesn't
    // do anything
    using elements_allocator = Allocator;
    using element_initializer = intrp2::ElementInitializer;
  };

  struct factory_creation {
    using factory_classes =
        tmpl::map<tmpl::pair<::Trigger, tmpl::list<Triggers::Always>>,
                  tmpl::pair<::Event, event_list>>;
  };
};

SPECTRE_TEST_CASE("Unit.ParallelAlgorithms.Interpolation.Runtime.Initialize",
                  "[Unit]") {
  using component = mock_interpolation_target<Metavariables>;
  register_factory_classes_with_charm<Metavariables>();

  using TriggerAndEvents = EventsAndTriggers::TriggerAndEvents;
  std::vector<TriggerAndEvents> trigger_and_events_storage{};
  tmpl::for_each<event_list>([&trigger_and_events_storage](auto event_v) {
    using DerivedEvent = tmpl::type_from<decltype(event_v)>;
    std::vector<std::unique_ptr<::Event>> events{};
    events.emplace_back(std::make_unique<DerivedEvent>());
    TriggerAndEvents storage{std::make_unique<Triggers::Always>(),
                             std::move(events)};
    trigger_and_events_storage.emplace_back(std::move(storage));
  });

  ActionTesting::MockRuntimeSystem<Metavariables> runner{
      {EventsAndTriggers{std::move(trigger_and_events_storage)}}};
  ActionTesting::set_phase(make_not_null(&runner),
                           Parallel::Phase::Initialization);
  for (size_t i = 0; i < num_elements; i++) {
    const std::string array_index = "TestEvent" + std::to_string(i);
    ActionTesting::emplace_array_component<component>(
        make_not_null(&runner), ActionTesting::NodeId{0},
        ActionTesting::LocalCoreId{0}, array_index);
    auto& box = ActionTesting::get_databox<component>(make_not_null(&runner),
                                                      array_index);

    // Shouldn't have been initialized yet
    db::mutate<intrp2::Tags::DbAccess>(
        [](const gsl::not_null<std::unique_ptr<db::Access>*> access) {
          CHECK_FALSE(*access);
        },
        make_not_null(&box));

    ActionTesting::next_action<component>(make_not_null(&runner), array_index);

    // Should now be initialized
    db::mutate<intrp2::Tags::DbAccess>(
        [](const gsl::not_null<std::unique_ptr<db::Access>*> access) {
          CHECK(*access);

          // ** because access is a pointer to a unique_ptr
          CHECK(db::get<AccessTag>(**access));

          // For the next test
          access->reset();

          CHECK_FALSE(*access);
        },
        make_not_null(&box));

    tmpl::for_each<event_list>([&](auto event_v) {
      using DerivedEvent = tmpl::type_from<decltype(event_v)>;
      DerivedEvent event{};
      // Only testing the current array index
      if (event.name() != array_index) {
        return;
      }

      ActionTesting::simple_action<
          component, intrp2::Actions::InitializeInterpolationTarget>(
          make_not_null(&runner), array_index,
          std::make_unique<DerivedEvent>());
    });

    // Again, should now be initialized
    db::mutate<intrp2::Tags::DbAccess>(
        [](const gsl::not_null<std::unique_ptr<db::Access>*> access) {
          CHECK(*access);

          // ** because access is a pointer to a unique_ptr
          CHECK(db::get<AccessTag>(**access));
        },
        make_not_null(&box));
  }
}
}  // namespace
