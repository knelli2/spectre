// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <optional>
#include <vector>

#include "Helpers/ParallelAlgorithms/NonlinearSolver/Algorithm.hpp"
#include "Options/String.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/CharmMain.tpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "ParallelAlgorithms/NonlinearSolver/NewtonRaphson/NewtonRaphson.hpp"
#include "Utilities/TMPL.hpp"

namespace PUP {
class er;
}  // namespace PUP

namespace helpers = TestHelpers::NonlinearSolver;

namespace {

struct LinearSolverGroup {
  static std::string name() { return "LinearSolver"; }
  inline const static std::string help {"Options for the linear solver"};
};

struct NonlinearSolverGroup {
  static std::string name() { return "NewtonRaphson"; }
  inline const static std::string help {"Options for the nonlinear solver"};
};

template <typename OperandTag>
struct ApplyNonlinearOperator {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const int /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    db::mutate<::NonlinearSolver::Tags::OperatorAppliedTo<OperandTag>>(
        [](const auto Ax, const auto& x) { *Ax = cube(x) - x; },
        make_not_null(&box), get<OperandTag>(box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

template <typename OperandTag, typename FieldTag>
struct ApplyLinearizedOperator {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const int /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    db::mutate<LinearSolver::Tags::OperatorAppliedTo<OperandTag>>(
        [](const auto Ap, const auto& dx, const auto& x) {
          *Ap = (3. * square(x) - 1) * dx;
        },
        make_not_null(&box), get<OperandTag>(box), get<FieldTag>(box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct Metavariables {
  inline const static std::string help{
      "Test the Newton-Raphson nonlinear solver algorithm"};

  using nonlinear_solver = NonlinearSolver::newton_raphson::NewtonRaphson<
      Metavariables, helpers::fields_tag, NonlinearSolverGroup>;
  using linear_solver = LinearSolver::gmres::Gmres<
      Metavariables, typename nonlinear_solver::linear_solver_fields_tag,
      LinearSolverGroup, false,
      typename nonlinear_solver::linear_solver_source_tag>;

  template <typename OperandTag>
  using apply_nonlinear_operator = ApplyNonlinearOperator<OperandTag>;
  template <typename OperandTag, typename FieldTag>
  using apply_linearized_operator =
      ApplyLinearizedOperator<OperandTag, FieldTag>;

  using component_list = helpers::component_list<Metavariables>;
  using observed_reduction_data_tags =
      helpers::observed_reduction_data_tags<Metavariables>;
  static constexpr bool ignore_unrecognized_command_line_options = false;
  static constexpr auto default_phase_order = helpers::default_phase_order;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

}  // namespace

extern "C" void CkRegisterMainModule() {
  Parallel::charmxx::register_main_module<Metavariables>();
  Parallel::charmxx::register_init_node_and_proc({}, {});
}
