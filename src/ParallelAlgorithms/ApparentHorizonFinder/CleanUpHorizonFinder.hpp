// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <map>
#include <set>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace ah {
template <typename DbTags>
void clean_up_and_go_to_next_time(
    const gsl::not_null<db::DataBox<DbTags>*> box,
    const gsl::not_null<Storage::NumberAndId*> number_and_id) {
  const std::map<Storage::NumberAndId, Storage::VolumeDataAndCallbacks>&
      volume_data_and_callbacks = db::get<Tags::VolumeDataAndCallbacks>(*box);

  // Check if we have more times to do horizon finds at. If we do, get the next
  // time we have. Check if it's the next horizon find. If it is, set the new
  // number and id so we try to find that horizon
  std::optional<Storage::NumberAndId> new_number_and_id{};
  if (not volume_data_and_callbacks.empty()) {
    const auto& volume_storage = *volume_data_and_callbacks.begin();
    const Storage::NumberAndId& potential_new_number_and_id =
        volume_storage.first;
    const Storage::VolumeDataAndCallbacks& volume_data = volume_storage.second;

    if (number_and_id->number + 1 == potential_new_number_and_id.number) {
      new_number_and_id = potential_new_number_and_id;

      db::mutate<Tags::SurfaceAndPoints>(
          [&number_and_id, &new_number_and_id,
           &volume_data](const gsl::not_null<std::unordered_map<
                             Storage::NumberAndId, Storage::SurfaceAndPoints>*>
                             surface_and_points) {
            (*surface_and_points)[new_number_and_id.value()] =
                surface_and_points->at(*number_and_id);

            // TODO
            // QUESTION: How to handle frames here? It's not guaranteed that the
            // next horizon find is in the same frame as this one. Maybe do a
            // transformation of coordinates?
            surface_and_points->at(new_number_and_id.value()).frame =
                volume_data.frame;
          },
          box);
    }
  }

  // Remove volume vars, interpolated vars, current strahlkorper and points,
  // reset fast flow, and mark this number and id as completed
  db::mutate<Tags::InterpolatedVars, Tags::VolumeDataAndCallbacks,
             Tags::SurfaceAndPoints, Tags::FastFlow,
             Tags::CompletedTemporalIds>(
      [&number_and_id](
          const gsl::not_null<
              std::map<Storage::NumberAndId, Storage::InterpolatedVars>*>
              interpolated_vars,
          const gsl::not_null<
              std::map<Storage::NumberAndId, Storage::VolumeDataAndCallbacks>*>
              volume_storage,
          const gsl::not_null<std::unordered_map<Storage::NumberAndId,
                                                 Storage::SurfaceAndPoints>*>
              surface_and_points,
          const gsl::not_null<FastFlow*> fast_flow,
          const gsl::not_null<std::set<Storage::NumberAndId>*>
              local_completed_ids) {
        interpolated_vars->erase(*number_and_id);
        volume_storage->erase(*number_and_id);
        surface_and_points->erase(*number_and_id);
        fast_flow->reset_for_next_find();
        local_completed_ids->emplace(*number_and_id);

        // We're fairly confident that we won't be receiving volume data
        // for a horizon find 1000 horizon finds later
        while (local_completed_ids->size() > 1000) {
          local_completed_ids->erase(*std::prev(local_completed_ids->end()));
        }
      },
      box);

  // Set the new number and id if it exists
  *number_and_id = new_number_and_id.value_or(*number_and_id);
}
}  // namespace ah
