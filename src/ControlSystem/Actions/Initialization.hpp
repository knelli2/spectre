// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Initialization/InitialData.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

#include <ostream>
#include "Parallel/Printf.hpp"

namespace Initialization {
namespace Actions {
template <typename Metavariables, typename ControlSystem>
struct InitializeControlSystem {
  static constexpr size_t deriv_order = ControlSystem::deriv_order;

  using initialization_tags =
      tmpl::list<control_system::Tags::ControlSystemInputs<ControlSystem>>;
  // using initialization_tags =
  //    tmpl::list<::Tags::Averager<2>, ::Tags::Controller<2>,
  //               ::Tags::TimescaleTuner, ::Tags::PrintOutput>;

  // using initialization_tags_to_keep = initialization_tags;

  // this is defined in ControlTranslation. Move to Tags
  using tags_to_be_initialized =
      tmpl::list<control_system::Tags::Averager<deriv_order>,
                 control_system::Tags::Controller<deriv_order>,
                 control_system::Tags::TimescaleTuner,
                 control_system::Tags::PrintOutput>;
  using simple_tags =
      tmpl::append<tags_to_be_initialized, typename ControlSystem::simple_tags>;

  using compute_tags = tmpl::list<>;

  template <
      typename DbTagsList, typename... InboxTags,  // typename Metavariables,
      typename ArrayIndex, typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/, ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    // Move all the control system inputs into their own tags in the databox so
    // we can use them easily later
    auto control_sys_inputs =
        db::get_mutable_reference<
            control_system::Tags::ControlSystemInputs<ControlSystem>>(
            make_not_null(&box))
            .get_inputs();
    ::Initialization::mutate_assign<tags_to_be_initialized>(
        make_not_null(&box),
        std::move(tuples::get<control_system::Tags::Averager<deriv_order>>(
            control_sys_inputs)),
        std::move(tuples::get<control_system::Tags::Controller<deriv_order>>(
            control_sys_inputs)),
        std::move(tuples::get<control_system::Tags::TimescaleTuner>(
            control_sys_inputs)),
        std::move(tuples::get<control_system::Tags::PrintOutput>(
            control_sys_inputs)));

    // Assign the control system name
    db::mutate<control_system::Tags::ControlSystemName>(
        make_not_null(&box),
        [](const gsl::not_null<std::string*> tag0, const std::string& name) {
          *tag0 = name;
        },
        ControlSystem::name());

    // Set the initial time between triggers using the initial timescale
    const auto& tuner = db::get<control_system::Tags::TimescaleTuner>(box);
    const double current_timescale = tuner.current_timescale()[0];
    db::mutate<control_system::Tags::Controller<deriv_order>>(
        make_not_null(&box),
        [&current_timescale](
            const gsl::not_null<::Controller<deriv_order>*> controller) {
          controller->assign_time_between_triggers(current_timescale);
        });
    return std::forward_as_tuple(std::move(box));
  }

  template <
      typename DataBox, typename... InboxTags,  // typename Metavariables,
      typename ArrayIndex, typename ActionList, typename ParallelComponent,
      Requires<not tmpl::all<initialization_tags,
                             tmpl::bind<db::tag_is_retrievable, tmpl::_1,
                                        tmpl::pin<DataBox>>>::value> = nullptr>
  static std::tuple<DataBox&&> apply(
      DataBox& /*box*/, const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    ERROR(
        "Dependencies not fulfilled. Did you forget to terminate the phase "
        "after removing options?");
  }
};
}  // namespace Actions
}  // namespace Initialization
