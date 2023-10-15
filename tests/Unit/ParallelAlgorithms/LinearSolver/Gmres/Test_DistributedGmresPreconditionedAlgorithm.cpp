// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <vector>

#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/Interval.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Helpers/Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Helpers/ParallelAlgorithms/LinearSolver/DistributedLinearSolverAlgorithmTestHelpers.hpp"
#include "Helpers/ParallelAlgorithms/LinearSolver/LinearSolverAlgorithmTestHelpers.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/CharmMain.tpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "ParallelAlgorithms/LinearSolver/Richardson/Richardson.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace PUP {
class er;
}  // namespace PUP

namespace helpers = LinearSolverAlgorithmTestHelpers;
namespace helpers_distributed = DistributedLinearSolverAlgorithmTestHelpers;

namespace {

struct ParallelGmres {
  inline const static std::string help
      {"Options for the iterative linear solver"};
};

struct Preconditioner {
  inline const static std::string help {"Options for the preconditioner"};
};

struct Metavariables {
  static constexpr const char* const help{
      "Test the preconditioned GMRES linear solver algorithm on multiple "
      "elements"};
  static constexpr size_t volume_dim = 1;
  using system =
      TestHelpers::domain::BoundaryConditions::SystemWithoutBoundaryConditions<
          volume_dim>;

  using linear_solver =
      LinearSolver::gmres::Gmres<Metavariables, helpers_distributed::fields_tag,
                                 ParallelGmres, true>;
  using preconditioner = LinearSolver::Richardson::Richardson<
      typename linear_solver::operand_tag, Preconditioner,
      typename linear_solver::preconditioner_source_tag>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DomainCreator<1>, tmpl::list<domain::creators::Interval>>>;
  };

  static constexpr auto default_phase_order = helpers::default_phase_order;
  using component_list = helpers_distributed::component_list<Metavariables>;
  using observed_reduction_data_tags =
      helpers::observed_reduction_data_tags<Metavariables>;
  static constexpr bool ignore_unrecognized_command_line_options = false;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

}  // namespace

extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<Metavariables>();
  Parallel::charmxx::register_init_node_and_proc(
      {&domain::creators::register_derived_with_charm}, {});
}
