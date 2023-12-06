// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
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
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits.hpp"
#include "Utilities/TypeTraits/CreateGetTypeAliasOrDefault.hpp"
#include "Utilities/TypeTraits/CreateIsCallable.hpp"

/// \cond

namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace intrp::Events {
struct MarkAsInterpolation;
}  // namespace intrp::Events
/// \endcond

/// Holds Actions for Interpolator and InterpolationTarget.
namespace intrp::Actions {

template <typename Metavariables>
struct InitializeInterpolationTarget2 {
 private:
  using factory_classes =
      typename std::decay_t<Metavariables>::factory_creation::factory_classes;
  using all_events = tmpl::at<factory_classes, ::Event>;
  using interpolation_events = tmpl::filter<
      all_events,
      std::is_base_of<intrp::Events::MarkAsInterpolation, tmpl::_1>>;
  using storage_tags =
      db::wrap_tags_in<intrp::Tags::TargetBoxStorage, interpolation_events>;

 public:
  using simple_tags = tmpl::push_back<storage_tags, intrp::Tags::DbAccesses>;
  using compute_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const EventsAndTriggers& events_and_triggers =
        Parallel::get<::Tags::EventsAndTriggers>(cache);

    const auto do_something = [&box](const ::Event& event) {
      // The array component element that we are on is for a specific event, so
      // only initialize the DB Access for that event box
      if (event.name() != array_index) {
        return;
      }

      // But we need the compile-time type of the event to initialize things
      // properly, so we have to loop over all the compile-time events and down
      // cast to check the event type
      tmpl::for_each<interpolation_events>([&event,
                                            &box](auto for_each_event_v) {
        using EventType = tmpl::type_from<decltype(for_each_event_v)>;

        const EventType* const derived = dynamic_cast<const EventType*>(&event);
        if (derived != nullptr) {
          db::mutate<intrp::Tags::DbAccesses>(
              [](const gsl::not_null<
                     std::unordered_map<std::string, db::Access* const>*>
                     accesses,
                 const std::unordered_map<
                     std::string, db::compte_databox_type<
                                      intrp::Tags::target_db_tags<EventType>>>&
                     storage_boxes) {
                (*accesses)[array_index] =
                    &db::as_access(storage_boxes.at(array_index));
              },
              make_not_null(&box),
              db::get<intrp::Tags::TargetBoxStorage<EventType>>(box));
        }
      });
    };

    events_and_triggers.for_each_event(do_something);

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

// The purpose of the metafunctions in this namespace is to allow
// InterpolationTarget::compute_target_points to omit an initialize
// function and a compute_tags and simple_tags type alias if it
// doesn't add anything to the DataBox.
namespace initialize_interpolation_target_detail {

CREATE_GET_TYPE_ALIAS_OR_DEFAULT(compute_tags)
CREATE_GET_TYPE_ALIAS_OR_DEFAULT(simple_tags)
CREATE_IS_CALLABLE(initialize)
CREATE_IS_CALLABLE_V(initialize)

}  // namespace initialize_interpolation_target_detail

/// \ingroup ActionsGroup
/// \brief Initializes an InterpolationTarget
///
/// Uses: nothing
///
/// DataBox changes:
/// - Adds:
///   - `Tags::IndicesOfFilledInterpPoints<TemporalId>`
///   - `Tags::IndicesOfInvalidInterpPoints<TemporalId>`
///   - `Tags::PendingTemporalIds<TemporalId>`
///   - `Tags::TemporalIds<TemporalId>` if target is non-sequential
///     `Tags::CurrentTemporalId<TemporalId>` if target is sequential
///   - `Tags::CompletedTemporalIds<TemporalId>`
///   - `Tags::InterpolatedVars<InterpolationTargetTag,TemporalId>`
///   - `::Tags::Variables<typename
///                   InterpolationTargetTag::vars_to_interpolate_to_target>`
/// - Removes: nothing
/// - Modifies: nothing
///
/// For requirements on InterpolationTargetTag, see InterpolationTarget
template <typename Metavariables, typename InterpolationTargetTag>
struct InitializeInterpolationTarget {
  using is_sequential =
      typename InterpolationTargetTag::compute_target_points::is_sequential;
  using TemporalId = typename InterpolationTargetTag::temporal_id::type;
  using return_tag_list_initial = tmpl::list<
      Tags::IndicesOfFilledInterpPoints<TemporalId>,
      Tags::IndicesOfInvalidInterpPoints<TemporalId>,
      Tags::PendingTemporalIds<TemporalId>,
      tmpl::conditional_t<is_sequential::value,
                          Tags::CurrentTemporalId<TemporalId>,
                          Tags::TemporalIds<TemporalId>>,
      Tags::CompletedTemporalIds<TemporalId>,
      Tags::InterpolatedVars<InterpolationTargetTag, TemporalId>,
      ::Tags::Variables<
          typename InterpolationTargetTag::vars_to_interpolate_to_target>>;

  using simple_tags = tmpl::append<
      return_tag_list_initial,
      initialize_interpolation_target_detail::get_simple_tags_or_default_t<
          typename InterpolationTargetTag::compute_target_points,
          tmpl::list<>>>;
  using compute_tags = tmpl::append<
      initialize_interpolation_target_detail::get_compute_tags_or_default_t<
          typename InterpolationTargetTag::compute_target_points, tmpl::list<>>,
      typename InterpolationTargetTag::compute_items_on_target>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    if constexpr (
        initialize_interpolation_target_detail::is_initialize_callable_v<
            typename InterpolationTargetTag::compute_target_points,
            const gsl::not_null<db::DataBox<DbTagsList>*>,
            const Parallel::GlobalCache<Metavariables>&>) {
      InterpolationTargetTag::compute_target_points::initialize(
          make_not_null(&box), cache);
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

}  // namespace intrp::Actions
