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
  using initialization_tags =
      tmpl::list<::Tags::Averager<2>, ::Tags::Controller<2>,
                 ::Tags::TimescaleTuner, ::Tags::PrintOutput>;

  using initialization_tags_to_keep = initialization_tags;

  // this is defined in ControlTranslation. Move to Tags
  using simple_tags =
      tmpl::flatten<tmpl::list<typename ControlSystem::simple_tags>>;

  using compute_tags = tmpl::list<>;

  template <
      typename DbTagsList, typename... InboxTags,  // typename Metavariables,
      typename ArrayIndex, typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/, ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    //::Initialization::mutate_assign<
    //    tmpl::list<control_system::Tags::ControlSystemName>>(
    //    make_not_null(&box), ControlSystem::name());
    db::mutate<control_system::Tags::ControlSystemName>(
        make_not_null(&box),
        [](const gsl::not_null<std::string*> tag0, const std::string& name) {
          *tag0 = name;
        },
        ControlSystem::name());
    const auto& tuner = db::get<::Tags::TimescaleTuner>(box);
    const double current_timescale = tuner.current_timescale()[0];
    db::mutate<::Tags::Controller<2>>(
        make_not_null(&box),
        [&current_timescale](const gsl::not_null<::Controller<2>*> controller) {
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
