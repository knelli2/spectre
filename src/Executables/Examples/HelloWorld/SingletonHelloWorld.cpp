// Distributed under the MIT License.
// See LICENSE.txt for details.

// [executable_example_includes]
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/CharmMain.tpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"

namespace PUP {
class er;
}  // namespace PUP
// [executable_example_includes]

// [executable_example_options]
namespace OptionTags {
struct Name {
  using type = std::string;
  inline const static std::string help{"A name"};
};
}  // namespace OptionTags

namespace Tags {
struct Name : db::SimpleTag {
  using type = std::string;
  using option_tags = tmpl::list<OptionTags::Name>;

  static constexpr bool pass_metavariables = false;
  static std::string create_from_options(const std::string& name) {
    return name;
  }
};
}  // namespace Tags
// [executable_example_options]

// [executable_example_action]
namespace Actions {
struct PrintMessage {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/) {
    Parallel::printf("Hello %s from process %d on node %d!\n",
                     Parallel::get<Tags::Name>(cache), sys::my_proc(),
                     sys::my_node());
  }
};
}  // namespace Actions
// [executable_example_action]

// [executable_example_singleton]
template <class Metavariables>
struct HelloWorld {
  using const_global_cache_tags = tmpl::list<Tags::Name>;
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Execute, tmpl::list<>>>;
  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;
  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache);
};

template <class Metavariables>
void HelloWorld<Metavariables>::execute_next_phase(
    const Parallel::Phase /* next_phase */,
    Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
  Parallel::simple_action<Actions::PrintMessage>(
      Parallel::get_parallel_component<HelloWorld>(
          *Parallel::local_branch(global_cache)));
}
// [executable_example_singleton]

// [executable_example_metavariables]
struct Metavars {
  using component_list = tmpl::list<HelloWorld<Metavars>>;

  inline const static std::string help{
      "Say hello from a singleton parallel component."};

  static constexpr std::array<Parallel::Phase, 3> default_phase_order{
      {Parallel::Phase::Initialization, Parallel::Phase::Execute,
       Parallel::Phase::Exit}};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
// [executable_example_metavariables]

// [executable_example_charm_init]
extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<Metavars>();
  Parallel::charmxx::register_init_node_and_proc({}, {});
}
// [executable_example_charm_init]
