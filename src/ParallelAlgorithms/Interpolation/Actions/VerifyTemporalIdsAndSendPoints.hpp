// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <sstream>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/SendPointsToInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"

namespace intrp::Actions {

template <typename InterpolationTargetTag>
struct VerifyTemporalIdsAndSendPoints;

namespace detail {
template <typename InterpolationTargetTag, typename ParallelComponent,
          typename DbTags, typename Metavariables>
void verify_temporal_ids_and_send_points_time_independent(
    const gsl::not_null<db::DataBox<DbTags>*> box,
    Parallel::GlobalCache<Metavariables>& cache) {
  using TemporalId = typename InterpolationTargetTag::temporal_id::type;
  std::stringstream ss{};
  ss << std::setprecision(std::numeric_limits<double>::digits10 + 4)
     << std::scientific;
  const ::Verbosity& verbosity = Parallel::get<intrp::Tags::Verbosity>(cache);
  const bool verbose_print = verbosity >= ::Verbosity::Verbose;
  ss << pretty_type::name<InterpolationTargetTag>() << ": ";

  // Move all PendingTemporalIds to TemporalIds, provided
  // that they are not already there, and fill new_temporal_ids
  // with the temporal_ids that were so moved.
  std::vector<TemporalId> new_temporal_ids{};
  db::mutate_apply<tmpl::list<Tags::TemporalIds<TemporalId>,
                              Tags::PendingTemporalIds<TemporalId>>,
                   tmpl::list<Tags::CompletedTemporalIds<TemporalId>>>(
      [&new_temporal_ids](
          const gsl::not_null<std::deque<TemporalId>*> ids,
          const gsl::not_null<std::deque<TemporalId>*> pending_ids,
          const std::deque<TemporalId>& completed_ids) {
        // This sort is needed because the ordering of these ids and pending ids
        // are not guaranteed. They aren't guaranteed because the elements must
        // send a communication to the this target with the the temporal id. So
        // it is possible that there are two interpolations that need to happen,
        // the earlier temporal id is sent first and the later temporal id is
        // sent second. But the later temporal id could arrive to this target
        // first because of charm communication latency. If this was it, then
        // there's nothing we could do and the later interpolation would happen
        // before the earlier one (even for sequential targets). However, there
        // is a scenario where this happens, except there is an interpolation in
        // progress so it doesn't send points to the interpolator when it
        // receives the later time first. In that case, we then receive the
        // earlier time second, and sort them here so they are in temporal
        // order. Note that this isn't a permanent actual fix, but so far works
        // in practice.
        alg::sort(*ids);
        alg::sort(*pending_ids);
        for (auto& id : *pending_ids) {
          if (std::find(completed_ids.begin(), completed_ids.end(), id) ==
                  completed_ids.end() and
              std::find(ids->begin(), ids->end(), id) == ids->end()) {
            ids->push_back(id);
            new_temporal_ids.push_back(id);
          }
        }
        pending_ids->clear();
      },
      box);
  if (InterpolationTargetTag::compute_target_points::is_sequential::value) {
    // Sequential: start interpolation only for the first new_temporal_id.
    if (not new_temporal_ids.empty()) {
      auto& my_proxy =
          Parallel::get_parallel_component<ParallelComponent>(cache);
      Parallel::simple_action<
          Actions::SendPointsToInterpolator<InterpolationTargetTag>>(
          my_proxy, new_temporal_ids.front());

      if (verbose_print) {
        ss << "Calling simple action to send points to interpolator at "
              "temporal id "
           << new_temporal_ids.front();
      }
    } else if (verbose_print) {
      ss << "No temporal ids to send points at.";
    }
  } else {
    // Non-sequential: start interpolation for all new_temporal_ids.
    auto& my_proxy = Parallel::get_parallel_component<ParallelComponent>(cache);
    for (const auto& id : new_temporal_ids) {
      Parallel::simple_action<
          Actions::SendPointsToInterpolator<InterpolationTargetTag>>(my_proxy,
                                                                     id);
    }
    if (verbose_print) {
      using ::operator<<;
      ss << "Calling simple action to send points to interpolator at temporal "
            "ids "
         << new_temporal_ids;
    }
  }

  if (verbose_print) {
    Parallel::printf("%s\n", ss.str());
  }
}

template <typename InterpolationTargetTag, typename ParallelComponent,
          typename DbTags, typename Metavariables, typename ArrayIndex>
void verify_temporal_ids_and_send_points_time_dependent(
    const gsl::not_null<db::DataBox<DbTags>*> box,
    Parallel::GlobalCache<Metavariables>& cache,
    const ArrayIndex& array_index) {
  using TemporalId = typename InterpolationTargetTag::temporal_id::type;
  std::stringstream ss{};
  ss << std::setprecision(std::numeric_limits<double>::digits10 + 4)
     << std::scientific;
  const ::Verbosity& verbosity = Parallel::get<intrp::Tags::Verbosity>(cache);
  const bool verbose_print = verbosity >= ::Verbosity::Verbose;
  ss << pretty_type::name<InterpolationTargetTag>() << ": ";

  const auto& pending_temporal_ids =
      db::get<Tags::PendingTemporalIds<TemporalId>>(*box);
  if (pending_temporal_ids.empty()) {
    if (verbose_print) {
      ss << "No pending temporal ids to send points at.";
      Parallel::printf("%s\n", ss.str());
    }
    return;  // Nothing to do if there are no pending temporal_ids.
  }

  auto& this_proxy = Parallel::get_parallel_component<ParallelComponent>(cache);
  double min_expiration_time = std::numeric_limits<double>::max();
  const Parallel::ArrayComponentId array_component_id =
      Parallel::make_array_component_id<ParallelComponent>(array_index);
  const bool at_least_one_pending_temporal_id_is_ready =
      ::Parallel::mutable_cache_item_is_ready<domain::Tags::FunctionsOfTime>(
          cache, array_component_id,
          [&this_proxy, &pending_temporal_ids, &min_expiration_time](
              const std::unordered_map<
                  std::string,
                  std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
                  functions_of_time) -> std::unique_ptr<Parallel::Callback> {
            min_expiration_time =
                std::min_element(functions_of_time.begin(),
                                 functions_of_time.end(),
                                 [](const auto& a, const auto& b) {
                                   return a.second->time_bounds()[1] <
                                          b.second->time_bounds()[1];
                                 })
                    ->second->time_bounds()[1];
            for (const auto& pending_id : pending_temporal_ids) {
              if (InterpolationTarget_detail::get_temporal_id_value(
                      pending_id) <= min_expiration_time) {
                // Success: at least one pending_temporal_id is ok.
                return std::unique_ptr<Parallel::Callback>{};
              }
            }
            // Failure: none of the pending_temporal_ids are ok.
            // Even though the GlobalCache docs say to only return a
            // PerformAlgorithmCallback, we return a SimpleActionCallback here
            // because it was already like this, and the effort to change it is
            // too great right now. This is alright because this code should
            // never be executed because the functions of time should always be
            // valid for times sent to the interpolation target
            return std::unique_ptr<Parallel::Callback>(
                new Parallel::SimpleActionCallback<
                    VerifyTemporalIdsAndSendPoints<InterpolationTargetTag>,
                    decltype(this_proxy)>(this_proxy));
          });

  if (not at_least_one_pending_temporal_id_is_ready) {
    // A callback has been set so that VerifyTemporalIdsAndSendPoints will
    // be called by the GlobalCache when domain::Tags::FunctionsOfTime is
    // updated.  So we can exit now.
    if (verbose_print) {
      ss << "At least one temporal id is not ready.";
      Parallel::printf("%s\n", ss.str());
    }
    return;
  }

  // Move up-to-date PendingTemporalIds to TemporalIds, provided
  // that they are not already there, and fill new_temporal_ids
  // with the temporal_ids that were so moved.
  std::vector<TemporalId> new_temporal_ids{};
  db::mutate_apply<tmpl::list<Tags::TemporalIds<TemporalId>,
                              Tags::PendingTemporalIds<TemporalId>>,
                   tmpl::list<Tags::CompletedTemporalIds<TemporalId>>>(
      [&min_expiration_time, &new_temporal_ids](
          const gsl::not_null<std::deque<TemporalId>*> ids,
          const gsl::not_null<std::deque<TemporalId>*> pending_ids,
          const std::deque<TemporalId>& completed_ids) {
        // This sort is needed because the ordering of these ids and pending ids
        // are not guaranteed. They aren't guaranteed because the elements must
        // send a communication to the this target with the the temporal id. So
        // it is possible that there are two interpolations that need to happen,
        // the earlier temporal id is sent first and the later temporal id is
        // sent second. But the later temporal id could arrive to this target
        // first because of charm communication latency. If this was it, then
        // there's nothing we could do and the later interpolation would happen
        // before the earlier one (even for sequential targets). However, there
        // is a scenario where this happens, except there is an interpolation in
        // progress so it doesn't send points to the interpolator when it
        // receives the later time first. In that case, we then receive the
        // earlier time second, and sort them here so they are in temporal
        // order. Note that this isn't a permanent actual fix, but so far works
        // in practice.
        alg::sort(*ids);
        alg::sort(*pending_ids);
        for (auto it = pending_ids->begin(); it != pending_ids->end();) {
          if (InterpolationTarget_detail::get_temporal_id_value(*it) <=
                  min_expiration_time and
              std::find(completed_ids.begin(), completed_ids.end(), *it) ==
                  completed_ids.end() and
              std::find(ids->begin(), ids->end(), *it) == ids->end()) {
            ids->push_back(*it);
            new_temporal_ids.push_back(*it);
            it = pending_ids->erase(it);
          } else {
            ++it;
          }
        }
      },
      box);

  if (InterpolationTargetTag::compute_target_points::is_sequential::value) {
    // Sequential: start interpolation only for the first new_temporal_id.
    if (not new_temporal_ids.empty()) {
      auto& my_proxy =
          Parallel::get_parallel_component<ParallelComponent>(cache);
      Parallel::simple_action<
          Actions::SendPointsToInterpolator<InterpolationTargetTag>>(
          my_proxy, new_temporal_ids.front());
      if (verbose_print) {
        ss << "Calling simple action to send points to interpolator at "
              "temporal id "
           << new_temporal_ids.front();
      }
    } else if (verbose_print) {
      ss << "No temporal ids to send points at.";
    }
  } else {
    // Non-sequential: start interpolation for all new_temporal_ids.
    auto& my_proxy = Parallel::get_parallel_component<ParallelComponent>(cache);
    for (const auto& id : new_temporal_ids) {
      Parallel::simple_action<
          Actions::SendPointsToInterpolator<InterpolationTargetTag>>(my_proxy,
                                                                     id);
    }
    // If there are still pending temporal_ids, call
    // VerifyTemporalIdsAndSendPoints again, so that those pending
    // temporal_ids can be waited for.
    if (not db::get<Tags::PendingTemporalIds<TemporalId>>(*box).empty()) {
      Parallel::simple_action<
          VerifyTemporalIdsAndSendPoints<InterpolationTargetTag>>(my_proxy);
    }

    if (verbose_print) {
      using ::operator<<;
      ss << "Calling simple action to send points to interpolator at temporal "
            "ids "
         << new_temporal_ids;
    }
  }

  if (verbose_print) {
    Parallel::printf("%s\n", ss.str());
  }
}
}  // namespace detail

/// \ingroup ActionsGroup
/// \brief Sends points to an Interpolator for verified temporal_ids.
///
/// VerifyTemporalIdsAndSendPoints is invoked on an InterpolationTarget.
///
/// In more detail, does the following:
/// - If any map is time-dependent:
///   - Moves verified PendingTemporalIds to TemporalIds, where
///     verified means that the FunctionsOfTime in the GlobalCache
///     are up-to-date for that TemporalId.  If no PendingTemporalIds are
///     moved, then VerifyTemporalIdsAndSendPoints sets itself as a
///     callback in the GlobalCache so that it is called again when the
///     FunctionsOfTime are mutated.
///   - If the InterpolationTarget is sequential, invokes
///     intrp::Actions::SendPointsToInterpolator for the first TemporalId.
///     (when interpolation is complete,
///      intrp::Actions::InterpolationTargetReceiveVars will begin interpolation
///     on the next TemporalId)
///   - If the InterpolationTarget is not sequential, invokes
///     intrp::Actions::SendPointsToInterpolator for all valid TemporalIds,
///     and then if PendingTemporalIds is non-empty it invokes itself.
///
/// - If all maps are time-independent:
///   - Moves all PendingTemporalIds to TemporalIds
///   - If the InterpolationTarget is sequential, invokes
///     intrp::Actions::SendPointsToInterpolator for the first TemporalId.
///     (when interpolation is complete,
///      intrp::Actions::InterpolationTargetReceiveVars will begin interpolation
///     on the next TemporalId)
///   - If the InterpolationTarget is not sequential, invokes
///     intrp::Actions::SendPointsToInterpolator for all TemporalIds.
///
/// Uses:
/// - DataBox:
///   - `intrp::Tags::PendingTeporalIds`
///   - `intrp::Tags::TeporalIds`
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies:
///   - `intrp::Tags::PendingTeporalIds`
///   - `intrp::Tags::TeporalIds`
///
/// For requirements on InterpolationTargetTag, see InterpolationTarget
template <typename InterpolationTargetTag>
struct VerifyTemporalIdsAndSendPoints {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index) {
    if constexpr (std::is_same_v<typename InterpolationTargetTag::
                                     compute_target_points::frame,
                                 ::Frame::Grid>) {
      detail::verify_temporal_ids_and_send_points_time_independent<
          InterpolationTargetTag, ParallelComponent>(make_not_null(&box),
                                                     cache);
    } else {
      const auto& domain =
          get<domain::Tags::Domain<Metavariables::volume_dim>>(cache);
      if (domain.is_time_dependent()) {
        if constexpr (Parallel::is_in_mutable_global_cache<
                          Metavariables, domain::Tags::FunctionsOfTime>) {
          detail::verify_temporal_ids_and_send_points_time_dependent<
              InterpolationTargetTag, ParallelComponent>(make_not_null(&box),
                                                         cache, array_index);
        } else {
          // We error here because the maps are time-dependent, yet
          // the cache does not contain FunctionsOfTime.  It would be
          // nice to make this a compile-time error; however, we want
          // the code to compile for the completely time-independent
          // case where there are no FunctionsOfTime in the cache at
          // all.  Unfortunately, checking whether the maps are
          // time-dependent is currently not constexpr.
          ERROR(
              "There is a time-dependent CoordinateMap in at least one "
              "of the Blocks, but FunctionsOfTime are not in the "
              "GlobalCache.  If you intend to use a time-dependent "
              "CoordinateMap, please add FunctionsOfTime to the GlobalCache.");
        }
      } else {
        detail::verify_temporal_ids_and_send_points_time_independent<
            InterpolationTargetTag, ParallelComponent>(make_not_null(&box),
                                                       cache);
      }
    }
  }
};
}  // namespace intrp::Actions
