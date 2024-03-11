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
#include "Evolution/EventsAndDenseTriggers/EventsAndDenseTriggers.hpp"
#include "Evolution/EventsAndDenseTriggers/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Events/MarkAsInterpolation.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Metavariables.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits.hpp"
#include "Utilities/TypeTraits/CreateGetTypeAliasOrDefault.hpp"
#include "Utilities/TypeTraits/CreateIsCallable.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
// IWYU pragma: no_forward_declare db::DataBox
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace Tags {
struct EventsAndTriggers;
}  // namespace Tags
/// \endcond

namespace intrp2 {

template <typename InterpolationEvents, typename BoxType>
void initialize_element_from_events(const gsl::not_null<BoxType*> box,
                                    const std::string& array_index,
                                    const ::Event& event) {
  static_assert(tt::is_a_v<tmpl::list, InterpolationEvents>);

  // We need the compile-time type of the event to initialize things properly,
  // so we have to loop over all the compile-time events and down cast to check
  // the event type
  tmpl::for_each<InterpolationEvents>([&](auto for_each_event_v) {
    using DerivedEventType = tmpl::type_from<decltype(for_each_event_v)>;
    const std::string& name = DerivedEventType{}.name();

    // The array component element that we are on is for a specific event, so
    // only initialize the DB Access for that event box
    if (array_index != name) {
      return;
    }

    const DerivedEventType* const derived =
        dynamic_cast<const DerivedEventType*>(&event);
    if (derived != nullptr) {
      // Yes this is strange to get a pointer to a unique pointer, but that's
      // what is required
      db::mutate<intrp2::Tags::DbAccess>(
          [&](const gsl::not_null<std::unique_ptr<db::Access>*> access) {
            const intrp2::Events::MarkAsInterpolation* const
                event_for_initialization =
                    dynamic_cast<const intrp2::Events::MarkAsInterpolation*>(
                        &event);

            ASSERT(event_for_initialization != nullptr,
                   "Event for initialization is nullptr");

            (*access) =
                event_for_initialization->initialize_target_element_box();
          },
          box);
    }
  });
}

struct ElementInitializer {
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables>
  static void apply(db::DataBox<DbTagList>& box,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const std::string& array_index) {
    using factory_classes =
        typename std::decay_t<Metavariables>::factory_creation::factory_classes;
    using all_events = tmpl::at<factory_classes, ::Event>;
    using interpolation_events = tmpl::filter<
        all_events,
        std::is_base_of<tmpl::pin<intrp2::Events::MarkAsInterpolation>,
                        tmpl::_1>>;

    const EventsAndTriggers& events_and_triggers =
        Parallel::get<::Tags::EventsAndTriggers>(cache);

    const auto for_each_event = [&](const ::Event& event) {
      if (array_index != event.name()) {
        return;
      }

      initialize_element_from_events<interpolation_events>(make_not_null(&box),
                                                           array_index, event);
    };

    events_and_triggers.for_each_event(for_each_event);
  }
};

/// Holds Actions for Interpolator and InterpolationTarget.
namespace Actions {
/*!
 * \brief Both an iterable and simple action for initializing the `db::Access`
 * for an interpolation target
 *
 */
struct InitializeInterpolationTarget {
  using const_global_cache_tags = tmpl::list<>;
  using simple_tags = tmpl::list<intrp2::Tags::DbAccess>;
  using compute_tags = tmpl::list<>;

  // Iterable action
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const std::string& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    using IntrpMetavars = typename Metavariables::intrp;
    static_assert(tt::assert_conforms_to_v<IntrpMetavars,
                                           intrp2::protocols::Metavariables>);
    using ElementInitializer = typename IntrpMetavars::element_initializer;

    ElementInitializer::template apply<ParallelComponent>(box, cache,
                                                          array_index);

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }

  // Simple action
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables>
  static void apply(db::DataBox<DbTagList>& box,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
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
    // cast to check the event type. If this array index doesn't correspond to
    // the event that was passed in, nothing will happen
    initialize_element_from_events<interpolation_events>(make_not_null(&box),
                                                         array_index, *event);
  }
};
}  // namespace Actions
}  // namespace intrp2
