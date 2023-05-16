// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Initialization/Tags.hpp"
#include "Evolution/Systems/Cce/IsBoundaryElement.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "IO/Logging/Tags.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Rational.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Cce::Actions {
template <typename Metavariables>
struct InitializeCcmDenseOutput {
  using simple_tags_from_options = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<logging::Tags::Verbosity<Cce::OptionTags::Cce>>;
  using mutable_global_cache_tags = tmpl::list<>;
  using return_tags =
      tmpl::list<::Tags::Variables<tmpl::list<Cce::Tags::TemporaryBondiJ>>,
                 ::Tags::Variables<typename Metavariables::ccm_psi0>>;
  using argument_tags = tmpl::list<Cce::Tags::LMax>;
  using compute_tags = tmpl::list<>;
  using simple_tags =
      tmpl::list<::Tags::Variables<tmpl::list<Cce::Tags::TemporaryBondiJ>>,
                 Cce::Tags::ExpectedNumberOfBoundaryElements,
                 ::Tags::Variables<typename Metavariables::ccm_psi0>>;

  static void apply(
      const gsl::not_null<::Variables<tmpl::list<Cce::Tags::TemporaryBondiJ>>*>
          temp_bondi_j,
      const gsl::not_null<::Variables<typename Metavariables::ccm_psi0>*>
          psi0_vars,
      const size_t l_max) {
    const size_t number_of_angular_grid_points =
        Spectral::Swsh::number_of_swsh_collocation_points(l_max);
    temp_bondi_j->initialize(number_of_angular_grid_points, 0.0);
    psi0_vars->initialize(number_of_angular_grid_points, 0.0);
  }
};

struct AddToExpectedNumberOnCcm {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/) {
    db::mutate<Cce::Tags::ExpectedNumberOfBoundaryElements>(
        make_not_null(&box), [](const gsl::not_null<size_t*> expected_number) {
          *expected_number += 1;
        });
  }
};

template <typename CcmComponent>
struct RegisterBoundaryElementsWithCcm {
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    // Don't send anything if we aren't on the outer boundary element
    if (not Cce::is_outer_boundary_element<false>(make_not_null(&inboxes), box,
                                                  cache)) {
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    Parallel::simple_action<AddToExpectedNumberOnCcm>(
        Parallel::get_parallel_component<CcmComponent>(cache));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

template <typename Psi0TagList>
struct InitializeCcmTags;

template <typename... Psi0Tags>
struct InitializeCcmTags<tmpl::list<Psi0Tags...>> {
  using simple_tags_from_options = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<Cce::Tags::LMax>;
  using mutable_global_cache_tags = tmpl::list<>;
  using return_tags = tmpl::list<Psi0Tags...>;
  using argument_tags = tmpl::list<Cce::Tags::LMax>;
  using compute_tags = tmpl::list<>;
  using simple_tags = return_tags;

  static void apply(const gsl::not_null<typename Psi0Tags::type*>... psi0_vars,
                    const size_t l_max) {
    const size_t number_of_angular_grid_points =
        Spectral::Swsh::number_of_swsh_collocation_points(l_max);

    const auto initialize_psi0_vars = [&number_of_angular_grid_points](
                                          const auto& psi0_var) {
      get(*psi0_var).data().destructive_resize(number_of_angular_grid_points);
      get(*psi0_var).data() = 0.0;
    };

    EXPAND_PACK_LEFT_TO_RIGHT(initialize_psi0_vars(psi0_vars));
  }
};
}  // namespace Cce::Actions
