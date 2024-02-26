// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/AsAccess.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits.hpp"
#include "Utilities/TypeTraits/CreateGetTypeAliasOrDefault.hpp"
#include "Utilities/TypeTraits/CreateIsCallable.hpp"

/// \cond
// IWYU pragma: no_forward_declare db::DataBox
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace intrp2::Events {
struct MarkAsInterpolation;
}  // namespace intrp2::Events
/// \endcond

namespace intrp2 {

/// Holds Actions for Interpolator and InterpolationTarget.
namespace Actions {

struct InitializeInterpolationTarget {
  using simple_tags = tmpl::list<intrp2::Tags::DbAccesses>;
  using compute_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const EventsAndTriggers& events_and_triggers =
        Parallel::get<::Tags::EventsAndTriggers>(cache);

    using factory_classes =
        typename std::decay_t<Metavariables>::factory_creation::factory_classes;
    using all_events = tmpl::at<factory_classes, ::Event>;
    using interpolation_events = tmpl::filter<
        all_events,
        std::is_base_of<tmpl::pin<intrp2::Events::MarkAsInterpolation>,
                        tmpl::_1>>;

    const auto do_something = [&box](const ::Event& event) {
      // The array component element that we are on is for a specific event, so
      // only initialize the DB Access for that event box
      if (event.name() != array_index) {
        return;
      }

      // But we need the compile-time type of the event to initialize things
      // properly, so we have to loop over all the compile-time events and down
      // cast to check the event type
      tmpl::for_each<interpolation_events>([&event, &box, &array_index](
                                               auto for_each_event_v) {
        using EventType = tmpl::type_from<decltype(for_each_event_v)>;

        const EventType* const derived = dynamic_cast<const EventType*>(&event);
        if (derived != nullptr) {
          db::mutate<intrp2::Tags::DbAccesses>(
              [](const gsl::not_null<
                  std::unordered_map<std::string, std::unique_ptr<db::Access>>*>
                     accesses) {
                const intrp2::Events::MarkAsInterpolation* const
                    event_for_initialization = dynamic_cast<
                        const intrp2::Events::MarkAsInterpolation*>(&event);

                (*accesses)[array_index] =
                    event_for_initialization->initialize_target_element_box();
              },
              make_not_null(&box));
        }
      });
    };

    events_and_triggers.for_each_event(do_something);

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

// TODO: Will need another simple action that does basically the exact same
// thing, except the simple action will need to also take the events that will
// be created at runtime.
struct InitializeInterpolationTargetCallback {
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables>
  static void apply(db::DataBox<DbTagList>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const std::string& array_index,
                    const std::unique_ptr<::Event>& event) {
    using factory_classes =
        typename std::decay_t<Metavariables>::factory_creation::factory_classes;
    using all_events = tmpl::at<factory_classes, ::Event>;
    using interpolation_events = tmpl::filter<
        all_events,
        std::is_base_of<tmpl::pin<intrp2::Events::MarkAsInterpolation>,
                        tmpl::_1>>;

    // We need the compile-time type of the event to initialize things
    // properly, so we have to loop over all the compile-time events and down
    // cast to check the event type
    tmpl::for_each<interpolation_events>([&event, &box,
                                          &array_index](auto for_each_event_v) {
      using EventType = tmpl::type_from<decltype(for_each_event_v)>;

      const EventType* const derived =
          dynamic_cast<const EventType*>(event.get());
      if (derived != nullptr) {
        db::mutate<intrp2::Tags::DbAccesses>(
            [](const gsl::not_null<
                std::unordered_map<std::string, std::unique_ptr<db::Access>>*>
                   accesses) {
              const intrp2::Events::MarkAsInterpolation* const
                  event_for_initialization =
                      dynamic_cast<const intrp2::Events::MarkAsInterpolation*>(
                          event.get());

              (*accesses)[array_index] =
                  event_for_initialization->initialize_target_element_box();
            },
            make_not_null(&box));
      }
    });
  }
};
}  // namespace Actions
}  // namespace intrp2
