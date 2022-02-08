// Distributed under the MIT License.
// See LICENSE.txt for details.

#define CATCH_CONFIG_RUNNER

#include "Framework/TestingFramework.hpp"

#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/Algorithms/AlgorithmGroup.hpp"
#include "Parallel/Algorithms/AlgorithmNodegroup.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InboxInserters.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Main.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "Parallel/Printf.hpp"
#include "Parallel/Tags/ResourceInfo.hpp"
#include "Parallel/TypeTraits.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

static constexpr int number_of_1d_array_elements = 14;

namespace PUP {
class er;
}  // namespace PUP
namespace db {
template <typename TagsList>
class DataBox;
}  // namespace db

struct TestMetavariables;

template <class Metavariables>
struct SingletonParallelComponent;

template <class Metavariables>
struct SingletonParallelComponentExclusive;

template <class Metavariables>
struct ArrayParallelComponent;

template <class Metavariables>
struct GroupParallelComponent;

template <class Metavariables>
struct NodegroupParallelComponent;

namespace SingletonActions {
struct Initialize {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    static_assert(Parallel::is_chare_type<Parallel::Algorithms::Singleton,
                                          ParallelComponent>::value,
                  "The ParallelComponent is not deduced to be the right type");
    // This test is run on 4 procs.
    if constexpr (std::is_same_v<
                      ParallelComponent,
                      SingletonParallelComponentExclusive<Metavariables>>) {
      // The exclusive proc is the last one, proc 3
      // Parallel::printf("On exclusive singleton. My proc = %d\n",
      //                 sys::my_proc());
      SPECTRE_PARALLEL_REQUIRE(sys::my_proc() == 3);
    } else {
      // The regular singleton is placed on the first available proc, which is 1
      // because we are avoiding the 0th proc
      // Parallel::printf("On normal singleton. My proc = %d\n",
      // sys::my_proc());
      SPECTRE_PARALLEL_REQUIRE(sys::my_proc() == 1);
    }
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace SingletonActions

namespace ArrayActions {
struct Initialize {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    static_assert(std::is_same_v<ParallelComponent,
                                 ArrayParallelComponent<TestMetavariables>>,
                  "The ParallelComponent is not deduced to be the right type");
    // We are running on 4 procs. We are avoiding the 0th proc, and one
    // singleton has requested to reserve the last proc (3).
    const int my_proc = sys::my_proc();
    // Parallel::printf("On dg-element. My proc = %d\n", sys::my_proc());
    SPECTRE_PARALLEL_REQUIRE(my_proc != 0);
    SPECTRE_PARALLEL_REQUIRE(my_proc != 3);
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace ArrayActions

namespace GroupActions {
struct Initialize {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    static_assert(std::is_same_v<ParallelComponent,
                                 GroupParallelComponent<TestMetavariables>>,
                  "The ParallelComponent is not deduced to be the right type");
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace GroupActions

namespace NodegroupActions {
struct Initialize {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    static_assert(std::is_same_v<ParallelComponent,
                                 NodegroupParallelComponent<TestMetavariables>>,
                  "The ParallelComponent is not deduced to be the right type");
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace NodegroupActions

template <class Metavariables>
struct SingletonParallelComponent {
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<SingletonActions::Initialize,
                 Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;
  using initialization_tags = Parallel::get_initialization_tags<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const typename Metavariables::Phase /*next_phase*/,
      const Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}
};

// [reserve_singleton_parallel_component]
template <class Metavariables>
struct SingletonParallelComponentExclusive {
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<SingletonActions::Initialize,
                 Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;
  using initialization_tags =
      tmpl::append<Parallel::get_initialization_tags<
                       Parallel::get_initialization_actions_list<
                           phase_dependent_action_list>>,
                   tmpl::list<Parallel::Tags::SingletonInfo<
                       SingletonParallelComponentExclusive>>>;

  static void execute_next_phase(
      const typename Metavariables::Phase /*next_phase*/,
      const Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}
};
// [reserve_singleton_parallel_component]

template <class Metavariables>
struct ArrayParallelComponent {
  using chare_type = Parallel::Algorithms::Array;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<ArrayActions::Initialize,
                 Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;
  using initialization_tags =
      tmpl::append<Parallel::get_initialization_tags<
                       Parallel::get_initialization_actions_list<
                           phase_dependent_action_list>>,
                   tmpl::list<Parallel::Tags::AvoidProc0>>;
  using array_index = int;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<initialization_tags>&
      /*initialization_items*/,
      const std::unordered_set<size_t>& procs_to_ignore = {}) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    auto& array_proxy =
        Parallel::get_parallel_component<ArrayParallelComponent>(local_cache);

    const size_t number_of_procs = static_cast<size_t>(sys::number_of_procs());
    size_t which_proc = 0;
    for (int i = 0; i < number_of_1d_array_elements; ++i) {
      while (procs_to_ignore.find(which_proc) != procs_to_ignore.end()) {
        which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
      }
      array_proxy[i].insert(global_cache, {}, which_proc);
      which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
    }
    array_proxy.doneInserting();
  }

  static void execute_next_phase(
      const typename Metavariables::Phase /*next_phase*/,
      const Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}
};

template <class Metavariables>
struct GroupParallelComponent {
  using chare_type = Parallel::Algorithms::Group;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<GroupActions::Initialize,
                 Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;
  using initialization_tags = Parallel::get_initialization_tags<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const typename Metavariables::Phase /*next_phase*/,
      Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}
};

template <class Metavariables>
struct NodegroupParallelComponent {
  using chare_type = Parallel::Algorithms::Nodegroup;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<NodegroupActions::Initialize,
                 Initialization::Actions::RemoveOptionsAndTerminatePhase>>>;
  using initialization_tags = Parallel::get_initialization_tags<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const typename Metavariables::Phase /*next_phase*/,
      Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}
};

struct TestMetavariables {
  using component_list =
      tmpl::list<SingletonParallelComponent<TestMetavariables>,
                 SingletonParallelComponentExclusive<TestMetavariables>,
                 ArrayParallelComponent<TestMetavariables>,
                 GroupParallelComponent<TestMetavariables>,
                 NodegroupParallelComponent<TestMetavariables>>;

  static constexpr const char* const help{
      "Test core placement of the Algorithm"};
  static constexpr bool ignore_unrecognized_command_line_options = false;

  enum class Phase { Initialization, Exit };
  template <typename... Tags>
  static Phase determine_next_phase(
      const gsl::not_null<
          tuples::TaggedTuple<Tags...>*> /*phase_change_decision_data*/,
      const Phase& /*current_phase*/,
      const Parallel::CProxy_GlobalCache<TestMetavariables>& /*cache_proxy*/) {
    return Phase::Exit;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling, &setup_memory_allocation_failure_reporting,
    &disable_openblas_multithreading};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};

using charmxx_main_component = Parallel::Main<TestMetavariables>;

#include "Parallel/CharmMain.tpp"  // IWYU pragma: keep
