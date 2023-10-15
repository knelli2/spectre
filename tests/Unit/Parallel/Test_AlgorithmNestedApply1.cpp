// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"  // IWYU pragma: keep
#include "Options/String.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/CharmMain.tpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Main.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace PUP {
class er;
}  // namespace PUP
namespace db {
template <typename TagsList>
class DataBox;
}  // namespace db

struct error_call_single_action_from_action {
  template <typename ParallelComponent, typename... DbTags,
            typename Metavariables, typename ArrayIndex>
  static void apply(db::DataBox<tmpl::list<DbTags...>>& /*box*/,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/) {
    // [bad_recursive_call]
    auto& local_parallel_component = *Parallel::local(
        Parallel::get_parallel_component<ParallelComponent>(cache));
    Parallel::simple_action<error_call_single_action_from_action>(
        local_parallel_component);
    // [bad_recursive_call]
  }
};

template <class Metavariables>
struct Component {
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization, tmpl::list<>>>;
  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      const Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    if (next_phase == Parallel::Phase::Execute) {
      auto& local_cache = *Parallel::local_branch(global_cache);
      Parallel::simple_action<error_call_single_action_from_action>(
          *Parallel::local(
              Parallel::get_parallel_component<Component>(local_cache)));
    }
  }
};

struct TestMetavariables {
  using component_list = tmpl::list<Component<TestMetavariables>>;

  static constexpr std::array<Parallel::Phase, 3> default_phase_order{
      {Parallel::Phase::Initialization, Parallel::Phase::Execute,
       Parallel::Phase::Exit}};

  inline const static std::string help {"Executable for testing"};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<TestMetavariables>();
  Parallel::charmxx::register_init_node_and_proc({}, {});
}
