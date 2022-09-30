// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Options/Options.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Time/EvolutionOrdering.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace Tags {
struct Time;
struct TimeStepId;
}  // namespace Tags
namespace domain::Tags {
struct FunctionsOfTime;
}
/// \endcond

namespace evolution {
namespace Tags {
/*!
 * \ingroup EventsAndTriggersGroup
 * \brief Previous time at which the trigger activated.
 *
 * \details This tag is populated with the most recent time the current trigger
 * fired prior to the current activation.
 * \note This tag is only populated with a meaningful value for events invoked
 * from `EventsAndDenseTriggers`; during other actions, the tag should not be
 * used.
 */
struct PreviousTriggerTime : db::SimpleTag {
  using type = std::optional<double>;
};
}  // namespace Tags

/// \ingroup EventsAndTriggersGroup
/// Class that checks dense triggers and runs events
class EventsAndDenseTriggers {
 private:
  template <typename Event>
  struct get_tags {
    using type = typename Event::compute_tags_for_observation_box;
  };

 public:
  using ConstructionType =
      std::vector<std::pair<std::unique_ptr<DenseTrigger>,
                            std::vector<std::unique_ptr<Event>>>>;

 private:
  struct TriggerRecord {
    double next_check;
    std::optional<bool> is_triggered;
    size_t num_events_ready;
    std::unique_ptr<DenseTrigger> trigger;
    std::vector<std::unique_ptr<Event>> events;

    // NOLINTNEXTLINE(google-runtime-references)
    void pup(PUP::er& p);
  };
  using Storage = std::vector<TriggerRecord>;

 public:
  EventsAndDenseTriggers() = default;
  explicit EventsAndDenseTriggers(ConstructionType events_and_triggers);

  template <typename DbTags>
  double next_trigger(const db::DataBox<DbTags>& box);

  enum class TriggeringState { Ready, NeedsEvolvedVariables, NotReady };

  /// Check triggers fire and whether all events to run are ready.  If
  /// this function returns anything other than NotReady, then
  /// rerunning it will skip checks (except for
  /// Event::needs_evolved_variables) until a successful call to
  /// reschedule.
  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename ComponentPointer>
  TriggeringState is_ready(const db::DataBox<DbTags>& box,
                           Parallel::GlobalCache<Metavariables>& cache,
                           const ArrayIndex& array_index,
                           const ComponentPointer component);

  /// Run events associated with fired triggers.  This must be called
  /// after is_ready returns something other then NotReady.  Any
  /// repeated calls will be no-ops until a successful call to
  /// reschedule.
  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename ComponentPointer>
  void run_events(db::DataBox<DbTags>& box,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ArrayIndex& array_index,
                  const ComponentPointer component);

  /// Schedule the next check.  This must be called after run_events
  /// for the current check.  Returns `true` on success, `false` if
  /// insufficient data is available and the call should be retried
  /// later.
  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename ComponentPointer>
  bool reschedule(const db::DataBox<DbTags>& box,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ArrayIndex& array_index,
                  const ComponentPointer component);

  /// Add a new trigger and set of events.  This can only be called
  /// during initialization.
  void add_trigger_and_events(std::unique_ptr<DenseTrigger> trigger,
                              std::vector<std::unique_ptr<Event>> events);

  template <typename F>
  void for_each_event(F&& f) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

 private:
  template <typename DbTags>
  void initialize(const db::DataBox<DbTags>& box);

  bool initialized() const;

  Storage events_and_triggers_;
  double next_check_ = std::numeric_limits<double>::signaling_NaN();
  evolution_less<double> before_{};
};

template <typename DbTags>
double EventsAndDenseTriggers::next_trigger(const db::DataBox<DbTags>& box) {
  if (UNLIKELY(not initialized())) {
    initialize(box);
  }

  if (events_and_triggers_.empty()) {
    return before_.infinity();
  }

  return next_check_;
}

template <typename DbTags, typename Metavariables, typename ArrayIndex,
          typename ComponentPointer>
EventsAndDenseTriggers::TriggeringState EventsAndDenseTriggers::is_ready(
    const db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& cache,
    const ArrayIndex& array_index, const ComponentPointer component) {
  ASSERT(initialized(), "Not initialized");
  ASSERT(not events_and_triggers_.empty(),
         "Should not be calling is_ready with no triggers");
  auto debug_print = [&box, &array_index, &cache](const std::string& message) {
    if (not is_zeroth_element(array_index)) {
      return;
    }
    double min_expiration_time{std::numeric_limits<double>::max()};
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
    if (not functions_of_time.empty()) {
      for (const auto& [name, f_of_t] : functions_of_time) {
        if (f_of_t->time_bounds()[1] < min_expiration_time) {
          min_expiration_time = f_of_t->time_bounds()[1];
        }
      }
    }
    Parallel::printf(
        "EventsAndDenseTriggers is_ready: i=%s texp=%1.20f t=%1.20f "
        "tstep=%1.20f "
        "tsub=%1.20f "
        "dt=%1.20f: %s\n",
        array_index, min_expiration_time, db::get<::Tags::Time>(box),
        db::get<::Tags::TimeStepId>(box).step_time().value(),
        db::get<::Tags::TimeStepId>(box).substep_time().value(),
        db::get<::Tags::TimeStep>(box).value(), message);
  };
  debug_print("Start is_ready");

  for (auto& trigger_entry : events_and_triggers_) {
    if (trigger_entry.next_check != next_check_) {
      debug_print("trigger_entry.next_check != next_check_, continue");
      continue;
    }
    if (not trigger_entry.is_triggered.has_value()) {
      debug_print("not trigger_entry.is_triggered.has_value()");
      const auto is_triggered = trigger_entry.trigger->is_triggered(
          box, cache, array_index, component);
      if (not is_triggered.has_value()) {
        debug_print("not is_triggered.has_value(), return NotReady");
        return TriggeringState::NotReady;
      }

      trigger_entry.is_triggered = *is_triggered;
    }

    if (not *trigger_entry.is_triggered) {
      debug_print("not *trigger_entry.is_triggered, continue");
      continue;
    }

    debug_print("checking ready events");
    for (; trigger_entry.num_events_ready < trigger_entry.events.size();
         ++trigger_entry.num_events_ready) {
      if (not trigger_entry.events[trigger_entry.num_events_ready]->is_ready(
              box, cache, array_index, component)) {
        debug_print("event not ready, return NotReady");
        return TriggeringState::NotReady;
      }
    }
  }

  debug_print("loop over triggers");
  for (auto& trigger_entry : events_and_triggers_) {
    if (trigger_entry.is_triggered != std::optional{true}) {
      debug_print(
          "trigger_entry.is_triggered != std::optional{true}, continue");
      continue;
    }
    for (const auto& event : trigger_entry.events) {
      if (event->needs_evolved_variables()) {
        debug_print(
            "event->needs_evolved_variables(), return NeedsEvolvedVariables");
        return TriggeringState::NeedsEvolvedVariables;
      }
    }
  }

  return TriggeringState::Ready;
}

template <typename DbTags, typename Metavariables, typename ArrayIndex,
          typename ComponentPointer>
void EventsAndDenseTriggers::run_events(
    db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& cache,
    const ArrayIndex& array_index, const ComponentPointer component) {
  ASSERT(initialized(), "Not initialized");
  ASSERT(not events_and_triggers_.empty(),
         "Should not be calling run_events with no triggers");
  using compute_tags = tmpl::remove_duplicates<tmpl::filter<
      tmpl::flatten<tmpl::transform<
          tmpl::at<typename Metavariables::factory_creation::factory_classes,
                   Event>,
          get_tags<tmpl::_1>>>,
      db::is_compute_tag<tmpl::_1>>>;
  auto debug_print = [&box, &array_index, &cache](const std::string& message) {
    if (not is_zeroth_element(array_index)) {
      return;
    }
    double min_expiration_time{std::numeric_limits<double>::max()};
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
    if (not functions_of_time.empty()) {
      for (const auto& [name, f_of_t] : functions_of_time) {
        if (f_of_t->time_bounds()[1] < min_expiration_time) {
          min_expiration_time = f_of_t->time_bounds()[1];
        }
      }
    }
    Parallel::printf(
        "EventsAndDenseTriggers run_events: i=%s texp=%1.20f t=%1.20f "
        "tstep=%1.20f "
        "tsub=%1.20f "
        "dt=%1.20f: %s\n",
        array_index, min_expiration_time, db::get<::Tags::Time>(box),
        db::get<::Tags::TimeStepId>(box).step_time().value(),
        db::get<::Tags::TimeStepId>(box).substep_time().value(),
        db::get<::Tags::TimeStep>(box).value(), message);
  };

  debug_print("Starting run_events");
  for (auto& trigger_entry : events_and_triggers_) {
    if (trigger_entry.is_triggered == std::optional{true}) {
      db::mutate<::evolution::Tags::PreviousTriggerTime>(
          make_not_null(&box),
          [&trigger_entry](const gsl::not_null<std::optional<double>*>
                               previous_trigger_time) {
            *previous_trigger_time =
                trigger_entry.trigger->previous_trigger_time();
          });
      const auto observation_box = make_observation_box<compute_tags>(box);
      debug_print(
          "trigger_entry.is_triggered == std::optional{true}, running events");
      for (const auto& event : trigger_entry.events) {
        event->run(observation_box, cache, array_index, component);
      }
      db::mutate<::evolution::Tags::PreviousTriggerTime>(
          make_not_null(&box), [](const gsl::not_null<std::optional<double>*>
                                      previous_trigger_time) {
            *previous_trigger_time =
                std::numeric_limits<double>::signaling_NaN();
          });
    }
    // Outside above if statement
    debug_print("trigger_entry.is_triggered = false");
    trigger_entry.is_triggered = false;
  }
}

template <typename DbTags, typename Metavariables, typename ArrayIndex,
          typename ComponentPointer>
bool EventsAndDenseTriggers::reschedule(
    const db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& cache,
    const ArrayIndex& array_index, const ComponentPointer component) {
  ASSERT(initialized(), "Not initialized");
  ASSERT(not events_and_triggers_.empty(),
         "Should not be calling run_events with no triggers");
  auto debug_print = [&box, &array_index, &cache](const std::string& message) {
    if (not is_zeroth_element(array_index)) {
      return;
    }
    double min_expiration_time{std::numeric_limits<double>::max()};
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
    if (not functions_of_time.empty()) {
      for (const auto& [name, f_of_t] : functions_of_time) {
        if (f_of_t->time_bounds()[1] < min_expiration_time) {
          min_expiration_time = f_of_t->time_bounds()[1];
        }
      }
    }
    Parallel::printf(
        "EventsAndDenseTriggers reschedule: i=%s texp=%1.20f t=%1.20f "
        "tstep=%1.20f "
        "tsub=%1.20f "
        "dt=%1.20f: %s\n",
        array_index, min_expiration_time, db::get<::Tags::Time>(box),
        db::get<::Tags::TimeStepId>(box).step_time().value(),
        db::get<::Tags::TimeStepId>(box).substep_time().value(),
        db::get<::Tags::TimeStep>(box).value(), message);
  };

  debug_print("Start reschedule");
  double new_next_check = before_.infinity();
  for (auto& trigger_entry : events_and_triggers_) {
    if (trigger_entry.next_check == next_check_) {
      debug_print("trigger_entry.next_check == next_check_");
      const std::optional<double> next_check =
          trigger_entry.trigger->next_check_time(box, cache, array_index,
                                                 component);
      if (not next_check.has_value()) {
        debug_print("not next_check.has_value(), return false");
        return false;
      }
      trigger_entry.next_check = *next_check;
      trigger_entry.num_events_ready = 0;
    }
    // Outside above if statement
    debug_print("trigger_entry.is_triggered.reset()");
    trigger_entry.is_triggered.reset();
    if (before_(trigger_entry.next_check, new_next_check)) {
      debug_print("before_(trigger_entry.next_check, new_next_check)");
      new_next_check = trigger_entry.next_check;
    }
  }
  if (is_zeroth_element(array_index)) {
    Parallel::printf(
        "EventsAndDenseTriggers reschedule: i=%s t=%1.20f next_check_=%1.20f\n",
        array_index, db::get<::Tags::Time>(box), new_next_check);
  }

  next_check_ = new_next_check;
  return true;
}

template <typename F>
void EventsAndDenseTriggers::for_each_event(F&& f) const {
  for (const auto& trigger_and_events : events_and_triggers_) {
    for (const auto& event : trigger_and_events.events) {
      f(*event);
    }
  }
}

template <typename DbTags>
void EventsAndDenseTriggers::initialize(const db::DataBox<DbTags>& box) {
  before_ = evolution_less<double>{
      db::get<::Tags::TimeStepId>(box).time_runs_forward()};
  next_check_ = db::get<::Tags::Time>(box);
  for (auto& trigger_record : events_and_triggers_) {
    trigger_record.next_check = next_check_;
  }
}
}  // namespace evolution

template <>
struct Options::create_from_yaml<evolution::EventsAndDenseTriggers> {
  using type = evolution::EventsAndDenseTriggers;
  template <typename Metavariables>
  static type create(const Options::Option& options) {
    return type(
        options.parse_as<typename type::ConstructionType, Metavariables>());
  }
};
