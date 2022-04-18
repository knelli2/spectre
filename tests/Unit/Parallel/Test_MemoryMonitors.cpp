// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Matrix.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/Initialize.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/MemoryMonitor/Initialize.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/MemoryMonitor/Tags.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/MonitorMemory.hpp"
#include "ParallelAlgorithms/Events/MonitorMemory.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace {

struct TimeTag {
  using type = double;
};

template <typename Metavariables>
using memory_tag_list =
    typename Initialization::Actions::InitializeMemoryMonitor<
        Metavariables>::simple_tags;

template <typename Metavariables>
struct MockMemoryMonitor {
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = int;
  using component_being_mocked = MemoryMonitor<Metavariables>;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<
          ActionTesting::InitializeDataBox<memory_tag_list<Metavariables>>>>>;
};

template <typename Metavariables>
struct MockObserverWriter {
  using component_being_mocked = observers::ObserverWriter<Metavariables>;

  using const_global_cache_tags =
      tmpl::list<observers::Tags::ReductionFileName>;

  using initialize_action_list =
      tmpl::list<::Actions::SetupDataBox,
                 observers::Actions::InitializeWriter<Metavariables>>;
  using initialization_tags =
      Parallel::get_initialization_tags<initialize_action_list>;

  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockNodeGroupChare;
  using array_index = int;

  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        initialize_action_list>>;
};

template <typename Metavariables>
struct SingletonParallelComponent {
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = int;
  using metavariables = Metavariables;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};

template <typename Metavariables>
struct GroupParallelComponent {
  using chare_type = ActionTesting::MockGroupChare;
  using array_index = int;
  using metavariables = Metavariables;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};

template <typename Metavariables>
struct NodegroupParallelComponent {
  using chare_type = ActionTesting::MockNodeGroupChare;
  using array_index = int;
  using metavariables = Metavariables;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};

template <typename Metavariables>
struct ArrayParallelComponent {
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using metavariables = Metavariables;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};

// This component deserves special mention. It is supposed to be playing the
// role of the DgElementArray, the component where we run the MonitorMemory
// event. However, inside MonitorMemory, there is an `if constexpr` check if the
// component we are monitoring is an array. If we are monitoring an array, then
// Parallel::contribute_to_reduction is called. If we make this component a
// MockArray, the test fails to build because the ATF doesn't support
// reductions yet. To get around this, we make this component a MockSingleton
// because the event only needs to be run on one "element". We also remove it
// from the components to monitor. Once the ATF supports reductions, this can be
// changed to a MockArray.
template <typename Metavariables>
struct FakeDgElementArray {
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = int;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<
          tmpl::list<domain::Tags::Element<3>>>>>>;
};

struct TestMetavariables {
  using observed_reduction_data_tags = tmpl::list<>;

  using component_list =
      tmpl::list<MockMemoryMonitor<TestMetavariables>,
                 MockObserverWriter<TestMetavariables>,
                 SingletonParallelComponent<TestMetavariables>,
                 GroupParallelComponent<TestMetavariables>,
                 FakeDgElementArray<TestMetavariables>,
                 NodegroupParallelComponent<TestMetavariables>>;

  enum class Phase { Initialization, Monitor, Exit };

  void pup(PUP::er& /*p*/) {}
};

using metavars = TestMetavariables;
using mem_mon_comp = MockMemoryMonitor<metavars>;
using obs_writer_comp = MockObserverWriter<metavars>;
using sing_comp = SingletonParallelComponent<metavars>;
using group_comp = GroupParallelComponent<metavars>;
using nodegroup_comp = NodegroupParallelComponent<metavars>;
using dg_elem_comp = FakeDgElementArray<metavars>;

void test_tags() {
  using holder_tag =
      mem_monitor::Tags::MemoryHolder<MockMemoryMonitor<TestMetavariables>>;
  TestHelpers::db::test_simple_tag<holder_tag>("MemoryHolder");

  const std::string subpath =
      mem_monitor::subfile_name<GroupParallelComponent<TestMetavariables>>();

  CHECK(subpath == "/MemoryMonitors/GroupParallelComponent");
}

struct BadArrayChareMetavariables {
  using component_list = tmpl::list<ArrayParallelComponent<TestMetavariables>>;

  enum class Phase { Initialization, Monitor, Exit };
};

void test_event_construction() {
  CHECK_THROWS_WITH(
      ([]() {
        std::vector<std::string> misspelled_component{"GlabolCahce"};

        Events::MonitorMemory<1, TimeTag> event{misspelled_component,
                                                Options::Context{}, metavars{}};
      }()),
      Catch::Contains(
          "Cannot monitor memory usage of unknown parallel component"));

  CHECK_THROWS_WITH(
      ([]() {
        std::vector<std::string> array_component{"ArrayParallelComponent"};

        Events::MonitorMemory<2, TimeTag> event{
            array_component, Options::Context{}, BadArrayChareMetavariables{}};
      }()),
      Catch::Contains("Currently, the only Array parallel component allowed to "
                      "be monitored is the DgElementArray."));
}

void test_wrong_box() {
  CHECK_THROWS_WITH(
      ([]() {
        db::DataBox<tmpl::list<>> empty_box{};
        Parallel::GlobalCache<metavars> cache{};

        mem_monitor::ContributeMemoryData<group_comp>::template apply<
            mem_mon_comp>(empty_box, cache, 0, 0.0, 0_st, 0.0);
      }()),
      Catch::Contains("Expected the DataBox for the MemoryMonitor"));
}

template <typename Component>
void check_output(const std::string& filename, const size_t num_nodes) {
  h5::H5File<h5::AccessType::ReadOnly> read_file{filename + ".h5"};
  INFO("Checking output of " + pretty_type::name<Component>());
  const auto& dataset =
      read_file.get<h5::Dat>(mem_monitor::subfile_name<Component>());
  const std::vector<std::string>& legend = dataset.get_legend();

  size_t num_columns;
  if constexpr (Parallel::is_singleton_v<Component>) {
    // time, proc, size
    num_columns = 3;
  } else if constexpr (Parallel::is_group_v<Component>) {
    // time, size on node 0, size on node 1, ...,  max size, proc of max size,
    // avg per node
    num_columns = num_nodes + 4;
  } else {
    // time, size on node 0, size on node 1, ..., avg per node
    num_columns = num_nodes + 2;
  }

  CHECK(legend.size() == num_columns);

  const Matrix data = dataset.get_data();
  // We only wrote one line of data for all components
  CHECK(data.rows() == 1);
  CHECK(data.columns() == num_columns);
}  // namespace

template <typename Component>
void run_group_actions(ActionTesting::MockRuntimeSystem<metavars>& runner,
                       const std::string& filename, const int num_nodes,
                       const int num_procs, const double time) {
  INFO("Checking actions of " + pretty_type::name<Component>());
  size_t num_branches;
  if constexpr (Parallel::is_group_v<Component>) {
    num_branches = static_cast<size_t>(num_procs);
  } else {
    num_branches = static_cast<size_t>(num_nodes);
  }

  for (int i = 0; i < static_cast<int>(num_branches); i++) {
    ActionTesting::invoke_queued_simple_action<Component>(
        make_not_null(&runner), i);
  }

  CHECK(ActionTesting::number_of_queued_simple_actions<mem_mon_comp>(
            runner, 0) == num_branches);

  const auto& mem_holder_tag = ActionTesting::get_databox_tag<
      mem_mon_comp, mem_monitor::Tags::MemoryHolder<Component>>(runner, 0);

  // Before we invoke the actions, the tag should be empty because nothing has
  // added to it yet
  CHECK(mem_holder_tag.empty());

  for (size_t i = 0; i < num_branches; i++) {
    // unordered_map<time, unordered_map<node/proc, size>>
    // Skip 0, otherwise an access error will occur
    if (i != 0) {
      CHECK(mem_holder_tag.at(time).size() == i);
    }
    ActionTesting::invoke_queued_simple_action<mem_mon_comp>(
        make_not_null(&runner), 0);
  }

  // After we invoke the actions, the tag should be empty because the time
  // should have been erased
  CHECK(mem_holder_tag.empty());

  // The last action should have called a threaded action to write data
  CHECK(ActionTesting::number_of_queued_threaded_actions<obs_writer_comp>(
            runner, 0) == 1);
  ActionTesting::invoke_queued_threaded_action<obs_writer_comp>(
      make_not_null(&runner), 0);

  check_output<Component>(filename, static_cast<size_t>(num_nodes));
}

void test_monitor_memory_event() {
  const std::string outfile_name{"TestMemoryMonitorOutput"};
  // clean up just in case
  if (file_system::check_if_file_exists(outfile_name + ".h5")) {
    file_system::rm(outfile_name + ".h5", true);
  }

  // 4 mock nodes, 3 mock cores per node
  const size_t num_nodes = 4;
  const size_t num_procs_per_node = 3;
  const size_t num_procs = num_nodes * num_procs_per_node;
  ActionTesting::MockRuntimeSystem<metavars> runner{
      {outfile_name}, {}, std::vector<size_t>(num_nodes, num_procs_per_node)};

  // Place simple components first (ones without initialization stuff)
  ActionTesting::emplace_singleton_component<sing_comp>(
      make_not_null(&runner), ActionTesting::NodeId{1},
      ActionTesting::LocalCoreId{1});
  ActionTesting::emplace_group_component<group_comp>(make_not_null(&runner));
  ActionTesting::emplace_nodegroup_component<nodegroup_comp>(
      make_not_null(&runner));

  // FakeDgElementArray that's actually a singleton
  const ElementId<3> element_id{0};
  const Element<3> element{element_id, {}};
  ActionTesting::emplace_singleton_component_and_initialize<dg_elem_comp>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, {element});

  // ObserverWriter
  ActionTesting::emplace_nodegroup_component<obs_writer_comp>(
      make_not_null(&runner));
  ActionTesting::next_action<obs_writer_comp>(make_not_null(&runner), 0);
  ActionTesting::next_action<obs_writer_comp>(make_not_null(&runner), 0);

  // MemoryMonitor
  ActionTesting::emplace_singleton_component_and_initialize<mem_mon_comp>(
      make_not_null(&runner), ActionTesting::NodeId{1},
      ActionTesting::LocalCoreId{1}, {});

  runner.set_phase(TestMetavariables::Phase::Monitor);

  auto& cache = ActionTesting::cache<mem_mon_comp>(runner, 0);

  std::vector<std::string> components_to_monitor{};
  // Note that we don't monitor the global caches here. This is because the
  // entry methods of the global caches used to compute memory and send it to
  // the MemoryMonitor are not compatible with the ATF.
  tmpl::for_each<metavars::component_list>([&components_to_monitor](
                                               auto component_v) {
    using component = tmpl::type_from<decltype(component_v)>;
    // Can't test DG element array in this test
    if constexpr (not std::is_same_v<component, FakeDgElementArray<metavars>>) {
      components_to_monitor.emplace_back(pretty_type::name<component>());
    }
  });

  // Create event
  Events::MonitorMemory<3, TimeTag> monitor_memory{
      components_to_monitor, Options::Context{}, metavars{}};

  const auto& element_box =
      ActionTesting::get_databox<dg_elem_comp,
                                 tmpl::list<domain::Tags::Element<3>>>(runner,
                                                                       0);

  // Run the event. This will queue a lot of actions
  const double time = 1.4;
  monitor_memory(time, element_box, cache, 0,
                 std::add_pointer_t<dg_elem_comp>{});

  // Check how many simple actions are queued:
  // - MemoryMonitor: 1
  // - ObserverWriter: 4 (one per node)
  // - Singleton: 1
  // - Group: 12 (one per core)
  // - NodeGroup: 4 (one per node)
  // - FakeDgElementArray: 0 (Because we didn't add this to the vector of
  //                          components to monitor above)
  tmpl::for_each<metavars::component_list>([&runner](auto component_v) {
    using component = tmpl::type_from<decltype(component_v)>;
    if constexpr (not std::is_same_v<component, FakeDgElementArray<metavars>>) {
      if constexpr (Parallel::is_singleton_v<component>) {
        CHECK(ActionTesting::number_of_queued_simple_actions<component>(
                  runner, 0) == 1);
      } else if constexpr (Parallel::is_group_v<component>) {
        for (int i = 0; i < static_cast<int>(num_procs); i++) {
          CHECK(ActionTesting::number_of_queued_simple_actions<component>(
                    runner, i) == 1);
        }
      } else if constexpr (Parallel::is_nodegroup_v<component>) {
        for (int i = 0; i < static_cast<int>(num_nodes); i++) {
          CHECK(ActionTesting::number_of_queued_simple_actions<component>(
                    runner, i) == 1);
        }
      }
    } else {
      CHECK(ActionTesting::is_simple_action_queue_empty<component>(runner, 0));
    }
  });

  // First invoke the simple actions on the singletons.
  ActionTesting::invoke_queued_simple_action<mem_mon_comp>(
      make_not_null(&runner), 0);
  ActionTesting::invoke_queued_simple_action<sing_comp>(make_not_null(&runner),
                                                        0);

  // These immediately call the observer writer to write data to disk
  CHECK(ActionTesting::number_of_queued_threaded_actions<obs_writer_comp>(
            runner, 0) == 2);

  ActionTesting::invoke_queued_threaded_action<obs_writer_comp>(
      make_not_null(&runner), 0);
  ActionTesting::invoke_queued_threaded_action<obs_writer_comp>(
      make_not_null(&runner), 0);

  // Check that the data was written correctly
  check_output<mem_mon_comp>(outfile_name, num_nodes);
  check_output<sing_comp>(outfile_name, num_nodes);

  // Now for the groups and nodegroups
  using group_list = tmpl::list<obs_writer_comp, group_comp, nodegroup_comp>;
  tmpl::for_each<group_list>([&runner, &outfile_name, &num_nodes, &num_procs,
                              &time](auto component_v) {
    using component = tmpl::type_from<decltype(component_v)>;
    run_group_actions<component>(runner, outfile_name, num_nodes, num_procs,
                                 time);
  });

  // All actions should be completed now
  CHECK(ActionTesting::number_of_queued_simple_actions<mem_mon_comp>(runner,
                                                                     0) == 0);
  CHECK(ActionTesting::number_of_queued_threaded_actions<obs_writer_comp>(
            runner, 0) == 0);

  // clean up just in case
  if (file_system::check_if_file_exists(outfile_name + ".h5")) {
    file_system::rm(outfile_name + ".h5", true);
  }
}

SPECTRE_TEST_CASE("Unit.Parallel.MemoryMonitors", "[Unit][Parallel]") {
  test_tags();
  test_event_construction();
  test_wrong_box();
  test_monitor_memory_event();
}
}  // namespace
