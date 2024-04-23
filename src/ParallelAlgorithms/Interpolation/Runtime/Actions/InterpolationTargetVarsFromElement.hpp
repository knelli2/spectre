// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Metafunctions.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/BlockLogicalCoordinates.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::Tags {
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
/// \endcond

namespace intrp2 {
namespace Actions {
template <typename Target>
struct ReceiveVolumeTensors {
  using TemporalId = typename Target::temporal_id_tag::type;
  using all_target_tensors = metafunctions::all_tags_to_observe<Target>;

  template <typename ParallelComponent, typename DbTags, typename Metavariables>
  static void apply(
      db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& cache,
      const std::string& /*array_index*/, const TemporalId& temporal_id,
      const std::unordered_map<std::string, std::vector<DataVector>>&
          received_interpolated_vars,
      const std::vector<size_t>& global_offsets,
      const bool check_for_invalid_points,
      const bool vars_have_already_been_received = false) {
    constexpr size_t Dim = Metavariables::volume_dim;
    // Check if we already have completed interpolation at this
    // temporal_id.
    db::Access& access =
        db::get_mutable_reference<Tags::DbAccess>(make_not_null(&box));

    if (UNLIKELY(
            alg::found(db::get<Tags::CompletedTemporalIds<TemporalId>>(access),
                       temporal_id))) {
      // The code will get into this 'if' statement in the following
      // scenario:
      // - There is at least one interpolation point exactly on the
      //   boundary of two or more Elements, so that
      //   InterpolationTargetVarsFromElement is called more than once
      //   with data for the same interpolation point (this is ok,
      //   and add_received_variables handles this).
      // - The only Elements that have not yet called
      //   InterpolationTargetVarsFromElement for this temporal_id are
      //   those that have data only for duplicated interpolation
      //   points, and the InterpolationTarget has already received
      //   that data from other Elements.
      // In this case, the InterpolationTarget proceeds to do its
      // work because it has all the data it needs. There is now
      // one more condition needed for the scenario that gets
      // us inside this 'if':
      // - The InterpolationTarget has already completed its work at
      //   this temporal_id, and it has cleaned up its data structures
      //   for this temporal_id before all of the remaining calls to
      //   InterpolationTargetVarsFromElement have occurred at this
      //   temporal_id, and now we are in one of those remaining
      //   calls.
      //
      // If this scenario occurs, we just return. This is because the
      // InterpolationTarget is done and there is nothing left to do
      // at this temporal_id.  Note that if there were extra work to
      // do at this temporal_id, then CompletedTemporalIds would not
      // have an entry for this temporal_id.
      return;
    }

    const auto& current_ids =
        db::get<Tags::CurrentTemporalIds<TemporalId>>(access);
    const auto& points = db::get<Tags::Points<Dim>>(access);
    const std::string& frame = db::get<Tags::Frame>(access);
    const size_t number_of_points = get<0>(points).size();
    // Reference in case this is a NaN
    const double& invalid_points_fill_value =
        db::get<Tags::InvalidPointsFillValue>(access);

    // If this is the first communication of data we are receiving, set up all
    // the tensors and fill them with the invalid value. That way, once the
    // interpolation is finished we have already filled the invalid values
    if (current_ids.count(temporal_id) == 0) {
      db::mutate<Tags::CurrentTemporalIds<TemporalId>,
                 Tags::NumberOfFilledPoints<TemporalId>,
                 Tags::InterpolatedVars<TemporalId>>(
          [&](const gsl::not_null<std::unordered_set<TemporalId>*>
                  mutable_current_ids,
              const gsl::not_null<std::unordered_map<TemporalId, size_t>*>
                  number_of_filled_points,
              const gsl::not_null<std::unordered_map<
                  TemporalId,
                  std::unordered_map<std::string, std::vector<DataVector>>>*>
                  interpolated_vars) {
            mutable_current_ids->insert(temporal_id);
            (*number_of_filled_points)[temporal_id] = 0;

            // Fill invalid values for each component of the tensor
            // QUESTION: Disable FPEs??
            auto& this_temporal_id_vars = (*interpolated_vars)[temporal_id];
            for (const auto& [name, received_data] :
                 received_interpolated_vars) {
              const size_t num_components = received_data.size();
              auto& vector_of_components = this_temporal_id_vars[name];
              vector_of_components.resize(num_components);
              for (size_t i = 0; i < num_components; i++) {
                vector_of_components[i] =
                    DataVector{number_of_points, invalid_points_fill_value};
              }
            }
          },
          make_not_null(&access));
    }

    if (not vars_have_already_been_received) {
      db::mutate<Tags::NumberOfFilledPoints<TemporalId>,
                 Tags::InterpolatedVars<TemporalId>>(
          [&](const gsl::not_null<std::unordered_map<TemporalId, size_t>*>
                  number_of_filled_points,
              const gsl::not_null<std::unordered_map<
                  TemporalId,
                  std::unordered_map<std::string, std::vector<DataVector>>>*>
                  interpolated_vars) {
            number_of_filled_points->at(temporal_id) += global_offsets.size();

            std::unordered_map<std::string, std::vector<DataVector>>&
                this_temporal_id_vars = interpolated_vars->at(temporal_id);

            // For each tensor that we were sent
            for (const auto& [name, received_data] :
                 received_interpolated_vars) {
              std::vector<DataVector>& interpolated_data =
                  this_temporal_id_vars.at(name);
              const size_t num_components = received_data.size();

              ASSERT(interpolated_data.size() == num_components,
                     "For tensor " << name << ", received number of components "
                                   << num_components
                                   << " is not the same as the existing number "
                                      "of components "
                                   << interpolated_data.size());

              // For each component of the tensor that we were sent
              for (size_t i = 0; i < num_components; i++) {
                const DataVector& received_component_data = received_data[i];
                DataVector& component_data = interpolated_data[i];

                ASSERT(
                    received_component_data.size() == global_offsets.size(),
                    "For tensor "
                        << name << ", component " << i
                        << ", the number of points received "
                        << received_component_data.size()
                        << " is not the same as the number of global offsets "
                        << global_offsets.size());

                // For each point we were sent
                for (size_t j = 0; j < global_offsets.size(); j++) {
                  const size_t global_offset = global_offsets[j];
                  component_data[global_offset] = received_component_data[j];
                }
              }
            }
          },
          make_not_null(&access));
    }

    // Here we check which indices are invalid. This way we know how many points
    // to receive total.
    if (check_for_invalid_points) {
      db::mutate<Tags::InvalidIndices<TemporalId>>(
          [&points, &frame, &temporal_id,
           &cache](const gsl::not_null<
                   std::unordered_map<TemporalId, std::unordered_set<size_t>>*>
                       all_invalid_indices) {
            ASSERT(not all_invalid_indices.contains(temporal_id),
                   "Checking invalid indices of target "
                       << pretty_type::name<Target>() << " at time "
                       << temporal_id
                       << " twice. This should only happen once.");
            auto& invalid_indices = all_invalid_indices[temporal_id];

            const intrp2::BlockCoords<Dim> block_logical_coordinates =
                ::block_logical_coordinates_in_frame(cache, temporal_id, points,
                                                     frame);

            for (size_t i = 0; i < block_logical_coordinates.size(); ++i) {
              if (not block_logical_coordinates[i].has_value()) {
                invalid_indices.insert(i);
              }
            }
          },
          make_not_null(&access));
    }

    const auto& all_invalid_indices =
        db::get<Tags::InvalidIndices<TemporalId>>(access);
    size_t received_number_of_points =
        db::get<Tags::NumberOfFilledPoints<TemporalId>>(access).at(temporal_id);
    if (all_invalid_indices.contains(temporal_id)) {
      received_number_of_points += all_invalid_indices.at(temporal_id).size();
    }
    // We don't have all the points yet. Wait for more.
    if (received_number_of_points != number_of_points) {
      return;
    }

    // We have all the points we need. First set the proper references in the
    // box
    const std::unordered_map<std::string, std::vector<DataVector>>&
        interpolated_vars =
            db::get<Tags::InterpolatedVars<TemporalId>>(access).at(temporal_id);
    const size_t number_of_sets_of_points =
        db::get<Tags::NumberOfSetsOfPoints>(access);
    const size_t number_of_points_in_single_set =
        db::get<Tags::NumberOfFilledPoints<TemporalId>>(access) /
        number_of_sets_of_points;
    for (size_t i = 0; i < number_of_sets_of_points; i++) {
      tmpl::for_each<all_target_tensors>([&interpolated_vars,
                                          &access](auto tensor_tag_v) {
        using TensorTag = tmpl::type_from<decltype(tensor_tag_v)>;
        const std::string& tag_name = db::tag_name<TensorTag>();

        // Skip this tag if it's not in the vars we received
        if (not interpolated_vars.contains(tag_name)) {
          return;
        }

        // Set the tensors in the actual access to be non-owning. Have them
        // point to the data in the unordered map so we avoid a bunch of
        // copies
        const std::vector<DataVector>& all_tensor_components =
            interpolated_vars.at(tag_name);
        db::mutate<TensorTag>(
            [&all_tensor_components](
                const gsl::not_null<typename TensorTag::type*> box_tensor) {
              for (size_t i = 0; i < all_tensor_components.size(); i++) {
                // NOLINTBEGIN
                box_tensor->get(i).set_data_ref(
                    const_cast<DataVector&>(all_tensor_components[i]).data() +
                        i * number_of_points_in_single_set,
                    number_of_points_in_single_set);
                // NOLINTEND
              }
            },
            make_not_null(&access));
      });

      // Now call all the callbacks
      const auto& callbacks = db::get<Tags::Callbacks<Target>>(access);
      for (const auto& callback : callbacks) {
        callback->invoke(access, cache, temporal_id);
      }
    }

// Unset the references so that we don't have dangling pointers. We only do this
// in debug mode so we save time.
// QUESTION: Should we always do this? Never? Seems like bad form to have
// dangling pointers
#ifdef SPECTRE_DEBUG
    tmpl::for_each<all_target_tensors>(
        [&interpolated_vars, &access](auto tensor_tag_v) {
          using TensorTag = tmpl::type_from<decltype(tensor_tag_v)>;
          const std::string& tag_name = db::tag_name<TensorTag>();

          // Skip this tag if it's not in the vars we received
          if (not interpolated_vars.contains(tag_name)) {
            return;
          }

          // Unset all references
          db::mutate<TensorTag>(
              [](const gsl::not_null<typename TensorTag::type*> box_tensor) {
                for (size_t i = 0; i < box_tensor->size(); i++) {
                  box_tensor->get(i).clear();
                }
              },
              make_not_null(&access));
        });
#endif

    // We are now finished with the callbacks, so clean up everything
    db::mutate<Tags::CurrentTemporalIds<TemporalId>,
               Tags::NumberOfFilledPoints<TemporalId>,
               Tags::InvalidIndices<TemporalId>,
               Tags::CompletedTemporalIds<TemporalId>,
               Tags::InterpolatedVars<TemporalId>>(
        [&temporal_id](
            const gsl::not_null<std::unordered_set<TemporalId>*>
                mutable_current_ids,
            const gsl::not_null<std::unordered_map<TemporalId, size_t>*>
                number_of_filled_points,
            const gsl::not_null<
                std::unordered_map<TemporalId, std::unordered_set<size_t>>*>
                invalid_indices,
            const gsl::not_null<std::unordered_set<TemporalId>*> completed_ids,
            const gsl::not_null<std::unordered_map<
                TemporalId,
                std::unordered_map<std::string, std::vector<DataVector>>>*>
                interpolated_vars) {
          mutable_current_ids->erase(temporal_id);
          number_of_filled_points->erase(temporal_id);
          invalid_indices->erase(temporal_id);
          interpolated_vars->erase(temporal_id);
          completed_ids->emplace_front(temporal_id);
          while (completed_ids->size() > 1000) {
            completed_ids->pop_back();
          }
        },
        make_not_null(&access));
  }
};
}  // namespace Actions
}  // namespace intrp2
