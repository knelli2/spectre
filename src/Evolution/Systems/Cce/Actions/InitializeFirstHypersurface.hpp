// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/Cce/Components/WorldtubeBoundary.hpp"
#include "Evolution/Systems/Cce/Initialize/InitializeJ.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/ScriPlusValues.hpp"
#include "IO/Observer/Actions/GetLockPointer.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

#include "Parallel/Printf.hpp"

/// \cond
struct DebugToggle;
/// \endcond

namespace Cce {
namespace Actions {

/*!
 * \ingroup ActionsGroup
 * \brief Given initial boundary data for \f$J\f$ and \f$\partial_r J\f$,
 * computes the initial hypersurface quantities \f$J\f$ and gauge values.
 *
 * \details This action is to be called after boundary data has been received,
 * but before the time-stepping evolution loop. So, it should be either late in
 * an initialization phase or early (before a `Actions::Goto` loop or similar)
 * in the `Evolve` phase.
 *
 * Internally, this dispatches to the call function of
 * `Tags::InitializeJ`, which designates a hypersurface initial data generator
 * chosen by input file options, `InitializeGauge`, and
 * `InitializeScriPlusValue<Tags::InertialRetardedTime>` to perform the
 * computations. Refer to the documentation for those mutators for mathematical
 * details.
 *
 * \note This action accesses the base tag `Cce::Tags::InitializeJBase`,
 * trusting that a tag that inherits from that base tag is present in the box or
 * the global cache. Typically, this tag should be added by the worldtube
 * boundary component, as the type of initial data is decided by the type of the
 * worldtube boundary data.
 */
template <bool UsesPartiallyFlatCartesianCoordinates,
          typename BoundaryComponent>
struct InitializeFirstHypersurface {
  using const_global_cache_tags =
      tmpl::list<Tags::LMax, Tags::NumberOfRadialPoints>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    // In some contexts, this action may get re-run (e.g. self-start procedure)
    // In those cases, we do not want to alter the existing hypersurface data,
    // so we just exit. However, we do want to re-run the action each time
    // the self start 'reset's from the beginning
    if (Parallel::get<DebugToggle>(cache)) {
      Parallel::printf(
          "CCE InitializeFirstHypersurface: At time %f, starting hypersurface "
          "computation.\n",
          db::get<::Tags::TimeStepId>(box).substep_time());
    }

    if (db::get<::Tags::TimeStepId>(box).slab_number() > 0 or
        not db::get<::Tags::TimeStepId>(box).is_at_slab_boundary()) {
      if (Parallel::get<DebugToggle>(cache)) {
        Parallel::printf(
            "CCE InitializeFirstHypersurface: At time %f, not altering "
            "previous hypersurface computation. Continuing\n",
            db::get<::Tags::TimeStepId>(box).substep_time());
      }
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    if (Parallel::get<DebugToggle>(cache)) {
      Parallel::printf(
          "CCE InitializeFirstHypersurface: At time %f, actually doing "
          "hypersurface calculation.\n",
          db::get<::Tags::TimeStepId>(box).substep_time());
    }
    // some initialization schemes need the hdf5_lock so that they can read
    // their own input data from disk.
    auto hdf5_lock = Parallel::local_branch(
                         Parallel::get_parallel_component<
                             observers::ObserverWriter<Metavariables>>(cache))
                         ->template local_synchronous_action<
                             observers::Actions::GetLockPointer<
                                 observers::Tags::H5FileLock>>();
    if constexpr (tt::is_a_v<AnalyticWorldtubeBoundary, BoundaryComponent>) {
      db::mutate_apply<typename InitializeJ::InitializeJ<false>::mutate_tags,
                       typename InitializeJ::InitializeJ<false>::argument_tags>(
          db::get<Tags::InitializeJBase>(box), make_not_null(&box),
          make_not_null(hdf5_lock));
    } else {
      db::mutate_apply<
          typename InitializeJ::InitializeJ<
              UsesPartiallyFlatCartesianCoordinates>::mutate_tags,
          typename InitializeJ::InitializeJ<
              UsesPartiallyFlatCartesianCoordinates>::argument_tags>(
          db::get<Tags::InitializeJBase>(box), make_not_null(&box),
          make_not_null(hdf5_lock));
    }
    db::mutate_apply<InitializeScriPlusValue<Tags::InertialRetardedTime>>(
        make_not_null(&box), db::get<::Tags::TimeStepId>(box).substep_time());
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Actions
}  // namespace Cce
