// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>
#include <tuple>
#include <utility>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system {
namespace Actions {
/*!
 * \ingroup ControlSystemGroup
 * \ingroup InitializationGroup
 * \brief Initialize items related to the control system
 *
 * DataBox:
 * - Uses:
 *   - `control_system::Tags::ControlSystemInputs<ControlSystem>`
 * - Adds:
 *   - `control_system::Tags::Averager<deriv_order>`
 *   - `control_system::Tags::Controller<deriv_order>`
 *   - `control_system::Tags::TimescaleTuner`
 *   - `control_system::Tags::ControlError`
 *   - `control_system::Tags::ControlSystemName`
 *   - `control_system::Tags::WriteDataToDisk`
 * - Removes: Nothing
 * - Modifies:
 *   - `control_system::Tags::Controller<deriv_order>`
 *
 * \note This action relies on the `SetupDataBox` aggregated initialization
 * mechanism, so `Actions::SetupDataBox` must be present in the `Initialization`
 * phase action list prior to this action.
 */
template <typename Metavariables, typename ControlSystem>
struct Initialize {
  static constexpr size_t deriv_order = ControlSystem::deriv_order;

  using initialization_tags =
      tmpl::list<control_system::Tags::ControlSystemInputs<ControlSystem>,
                 control_system::Tags::WriteDataToDisk>;

  using initialization_tags_to_keep =
      tmpl::list<control_system::Tags::WriteDataToDisk>;

  using tags_to_be_initialized =
      tmpl::list<control_system::Tags::Averager<deriv_order - 1>,
                 control_system::Tags::Controller<deriv_order>,
                 control_system::Tags::TimescaleTuner,
                 control_system::Tags::ControlError<ControlSystem>,
                 control_system::Tags::ControlSystemName>;

  using simple_tags =
      tmpl::append<tags_to_be_initialized, typename ControlSystem::simple_tags>;

  using compute_tags = tmpl::list<>;

  using const_global_cache_tags =
      tmpl::list<control_system::Tags::RestrictToRotationAboutZAxis>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    // Move all the control system inputs into their own tags in the databox so
    // we can use them easily later
    const auto& option_holder =
        db::get<control_system::Tags::ControlSystemInputs<ControlSystem>>(box);
    ::Initialization::mutate_assign<tags_to_be_initialized>(
        make_not_null(&box), option_holder.averager, option_holder.controller,
        option_holder.tuner, option_holder.control_error,
        ControlSystem::name());

    // Set the initial time between updates and measurements
    const auto& measurement_timescales =
        get<control_system::Tags::MeasurementTimescales>(cache);
    const auto& measurement_timescale_func =
        *(measurement_timescales.at(ControlSystem::name()));
    const double initial_time = measurement_timescale_func.time_bounds()[0];
    const double measurement_timescale =
        min(measurement_timescale_func.func(initial_time)[0]);
    const auto& tuner = db::get<control_system::Tags::TimescaleTuner>(box);
    const double current_min_damp_timescale = min(tuner.current_timescale());
    db::mutate<control_system::Tags::Controller<deriv_order>,
               control_system::Tags::Averager<deriv_order - 1>>(
        make_not_null(&box),
        [&current_min_damp_timescale, &initial_time, &measurement_timescale](
            const gsl::not_null<::Controller<deriv_order>*> controller,
            const gsl::not_null<::Averager<deriv_order - 1>*> averager) {
          controller->set_initial_update_time(initial_time);
          controller->assign_time_between_updates(current_min_damp_timescale);
          averager->assign_time_between_measurements(measurement_timescale);
        });

    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace Actions
}  // namespace control_system
