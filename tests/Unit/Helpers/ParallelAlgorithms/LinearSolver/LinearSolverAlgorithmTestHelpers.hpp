// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DynamicMatrix.hpp"
#include "DataStructures/DynamicVector.hpp"
#include "IO/Observer/Actions/RegisterWithObservers.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "NumericalAlgorithms/Convergence/HasConverged.hpp"
#include "NumericalAlgorithms/Convergence/Reason.hpp"
#include "NumericalAlgorithms/Convergence/Tags.hpp"
#include "Options/String.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/Goto.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "ParallelAlgorithms/LinearSolver/Actions/MakeIdentityIfSkipped.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace LinearSolverAlgorithmTestHelpers {

namespace OptionTags {
struct LinearOperator {
  inline const static std::string help {"The linear operator A to invert."};
  using type = blaze::DynamicMatrix<double>;
};
struct Source {
  inline const static std::string help {"The source b in the equation Ax=b."};
  using type = blaze::DynamicVector<double>;
};
struct InitialGuess {
  inline const static std::string help {"The initial guess for the vector x."};
  using type = blaze::DynamicVector<double>;
};
struct ExpectedResult {
  inline const static std::string help {"The solution x in the equation Ax=b"};
  using type = blaze::DynamicVector<double>;
};
struct ExpectedConvergenceReason {
  static std::string name() { return "ConvergenceReason"; }
  inline const static std::string help {"The expected convergence reason"};
  using type = Convergence::Reason;
};
}  // namespace OptionTags

struct LinearOperator : db::SimpleTag {
  using type = blaze::DynamicMatrix<double>;
  using option_tags = tmpl::list<OptionTags::LinearOperator>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& linear_operator) {
    return linear_operator;
  }
};

struct Source : db::SimpleTag {
  using type = blaze::DynamicVector<double>;
  using option_tags = tmpl::list<OptionTags::Source>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& source) { return source; }
};

struct InitialGuess : db::SimpleTag {
  using type = blaze::DynamicVector<double>;
  using option_tags = tmpl::list<OptionTags::InitialGuess>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& initial_guess) {
    return initial_guess;
  }
};

struct ExpectedResult : db::SimpleTag {
  using type = blaze::DynamicVector<double>;
  using option_tags = tmpl::list<OptionTags::ExpectedResult>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& expected_result) {
    return expected_result;
  }
};

struct ExpectedConvergenceReason : db::SimpleTag {
  using type = Convergence::Reason;
  using option_tags = tmpl::list<OptionTags::ExpectedConvergenceReason>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option_value) {
    return option_value;
  }
};

// The vector `x` we want to solve for
struct VectorTag : db::SimpleTag {
  using type = blaze::DynamicVector<double>;
};

using fields_tag = VectorTag;

template <typename OperandTag>
struct ComputeOperatorAction {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const int /*array_index*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ActionList /*meta*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ParallelComponent* const /*meta*/) {
    db::mutate<LinearSolver::Tags::OperatorAppliedTo<OperandTag>>(
        [](const gsl::not_null<blaze::DynamicVector<double>*>
               operator_applied_to_operand,
           const blaze::DynamicMatrix<double>& linear_operator,
           const blaze::DynamicVector<double>& operand) {
          *operator_applied_to_operand = linear_operator * operand;
        },
        make_not_null(&box), get<LinearOperator>(box), get<OperandTag>(box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

// Checks for the correct solution after the algorithm has terminated.
template <typename OptionsGroup>
struct TestResult {
  using const_global_cache_tags =
      tmpl::list<ExpectedResult, ExpectedConvergenceReason>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const int /*array_index*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ActionList /*meta*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ParallelComponent* const /*meta*/) {
    const auto& has_converged =
        get<Convergence::Tags::HasConverged<OptionsGroup>>(box);
    SPECTRE_PARALLEL_REQUIRE(has_converged);
    SPECTRE_PARALLEL_REQUIRE(has_converged.reason() ==
                             get<ExpectedConvergenceReason>(box));
    const auto& result = get<VectorTag>(box);
    const auto& expected_result = get<ExpectedResult>(box);
    for (size_t i = 0; i < expected_result.size(); i++) {
      SPECTRE_PARALLEL_REQUIRE(result[i] == approx(expected_result[i]));
    }
    return {Parallel::AlgorithmExecution::Pause, std::nullopt};
  }
};

struct InitializeElement {
  using simple_tags = tmpl::list<VectorTag, ::Tags::FixedSource<VectorTag>>;
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const int /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    Initialization::mutate_assign<simple_tags>(
        make_not_null(&box), get<InitialGuess>(box), get<Source>(box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

namespace detail {

template <typename Preconditioner>
struct init_preconditioner_impl {
  using type = typename Preconditioner::initialize_element;
};

template <>
struct init_preconditioner_impl<void> {
  using type = tmpl::list<>;
};

template <typename Preconditioner>
using init_preconditioner =
    typename init_preconditioner_impl<Preconditioner>::type;

template <typename Preconditioner>
struct register_preconditioner_impl {
  using type = typename Preconditioner::register_element;
};

template <>
struct register_preconditioner_impl<void> {
  using type = tmpl::list<>;
};

template <typename Preconditioner>
using register_preconditioner =
    typename register_preconditioner_impl<Preconditioner>::type;

template <typename Preconditioner>
struct run_preconditioner_impl {
  using type =
      tmpl::list<ComputeOperatorAction<typename Preconditioner::fields_tag>,
                 // [make_identity_if_skipped]
                 typename Preconditioner::template solve<ComputeOperatorAction<
                     typename Preconditioner::operand_tag>>,
                 LinearSolver::Actions::MakeIdentityIfSkipped<Preconditioner>
                 // [make_identity_if_skipped]
                 >;
};

template <>
struct run_preconditioner_impl<void> {
  using type = tmpl::list<>;
};

template <typename Preconditioner>
using run_preconditioner =
    typename run_preconditioner_impl<Preconditioner>::type;

}  // namespace detail

template <typename Metavariables>
struct ElementArray {
  using chare_type = Parallel::Algorithms::Array;
  using array_index = int;
  using metavariables = Metavariables;
  using linear_solver = typename Metavariables::linear_solver;
  using preconditioner = typename Metavariables::preconditioner;

  // In each step of the algorithm we must provide A(p). The linear solver then
  // takes care of updating x and p, as well as the internal variables r, its
  // magnitude and the iteration step number.
  /// [action_list]
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<InitializeElement,
                     typename linear_solver::initialize_element,
                     ComputeOperatorAction<fields_tag>,
                     detail::init_preconditioner<preconditioner>,
                     Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Register,
          tmpl::list<typename linear_solver::register_element,
                     detail::register_preconditioner<preconditioner>,
                     Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Solve,
          tmpl::list<
              typename linear_solver::template solve<tmpl::list<
                  detail::run_preconditioner<preconditioner>,
                  ComputeOperatorAction<typename linear_solver::operand_tag>>>,
              Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<TestResult<typename linear_solver::options_group>>>>;
  /// [action_list]
  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;
  using const_global_cache_tags =
      tmpl::list<LinearOperator, Source, InitialGuess, ExpectedResult>;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<simple_tags_from_options>&
          initialization_items,
      const std::unordered_set<size_t>& /*procs_to_ignore*/ = {}) {
    auto& local_component = Parallel::get_parallel_component<ElementArray>(
        *Parallel::local_branch(global_cache));
    local_component[0].insert(global_cache, initialization_items, 0);
    local_component.doneInserting();
  }

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_component = Parallel::get_parallel_component<ElementArray>(
        *Parallel::local_branch(global_cache));
    local_component.start_phase(next_phase);
  }
};

// After the algorithm completes we perform a cleanup phase that checks the
// expected output file was written and deletes it.
template <bool CheckExpectedOutput, bool ExpectReductions, bool ExpectVolume>
struct CleanOutput {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& reductions_file_name =
        get<observers::Tags::ReductionFileName>(box) + ".h5";
    if (file_system::check_if_file_exists(reductions_file_name)) {
      file_system::rm(reductions_file_name, true);
    } else if (CheckExpectedOutput and ExpectReductions) {
      ERROR("Expected reductions file '" << reductions_file_name
                                         << "' does not exist");
    }
    const auto& volume_file_name =
        get<observers::Tags::VolumeFileName>(box) + "0.h5";
    if (file_system::check_if_file_exists(volume_file_name)) {
      file_system::rm(volume_file_name, true);
    } else if (CheckExpectedOutput and ExpectVolume) {
      ERROR("Expected volume file '" << volume_file_name << "' does not exist");
    }
    return {Parallel::AlgorithmExecution::Pause, std::nullopt};
  }
};

template <typename Metavariables, bool ExpectReductions = true,
          bool ExpectVolume = false>
struct OutputCleaner {
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::list<CleanOutput<false, ExpectReductions, ExpectVolume>>>,

      Parallel::PhaseActions<
          Parallel::Phase::Cleanup,
          tmpl::list<CleanOutput<true, ExpectReductions, ExpectVolume>>>>;
  using simple_tags_from_options = Parallel::get_simple_tags_from_options<
      Parallel::get_initialization_actions_list<phase_dependent_action_list>>;

  static void execute_next_phase(
      const Parallel::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_component = Parallel::get_parallel_component<OutputCleaner>(
        *Parallel::local_branch(global_cache));
    local_component.start_phase(next_phase);
  }
};

// [default_phase_order_array]
static constexpr std::array<Parallel::Phase, 6> default_phase_order{
    {Parallel::Phase::Initialization, Parallel::Phase::Register,
     Parallel::Phase::Solve, Parallel::Phase::Testing, Parallel::Phase::Cleanup,
     Parallel::Phase::Exit}};
// [default_phase_order_array]

namespace detail {

template <typename LinearSolver>
struct get_component_list_impl {
  using type = typename LinearSolver::component_list;
};

template <>
struct get_component_list_impl<void> {
  using type = tmpl::list<>;
};

template <typename LinearSolver>
using get_component_list = typename get_component_list_impl<LinearSolver>::type;

}  // namespace detail

template <typename Metavariables>
using component_list = tmpl::push_back<
    tmpl::append<
        detail::get_component_list<typename Metavariables::linear_solver>,
        detail::get_component_list<typename Metavariables::preconditioner>>,
    ElementArray<Metavariables>, observers::Observer<Metavariables>,
    observers::ObserverWriter<Metavariables>, OutputCleaner<Metavariables>>;

template <typename Metavariables>
using observed_reduction_data_tags =
    observers::collect_reduction_data_tags<tmpl::flatten<tmpl::list<
        typename Metavariables::linear_solver,
        tmpl::conditional_t<
            std::is_same_v<typename Metavariables::preconditioner, void>,
            tmpl::list<>, typename Metavariables::preconditioner>>>>;

}  // namespace LinearSolverAlgorithmTestHelpers
