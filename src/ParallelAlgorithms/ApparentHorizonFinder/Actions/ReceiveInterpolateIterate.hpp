// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/Callback.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Actions/InterpolateToTargetPoints.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/CleanUpHorizonFinder.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeTargetPoints.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Destination.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/SingleHorizonFindIteration.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Storage.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace ah {
namespace Actions {
// Important concepts about the two different communication/parallelization
// algorithms for horizon finding. Listing both pros and cons
//
// Interpolator/target method
//  - Every element on a core only communicates volume data to the interpolator
//    on their core. This avoids network communication of large buffer of data.
//  - Volume data is spread out over all cores. No worry about running out of
//    memory.
//  - Interpolation to target points is done in parallel on the interpolator
//    cores.
//  - Every element communicates to the target (single core) telling the target
//    to send points to the interpolator.
//  - Target sends same communication of points to all cores of the
//    interpolator.
//  - Every iteration of the horizon find is a communication from the target to
//    all interpolator cores, then a parallel interpolation onto target points,
//    then a communication from a subset of the interpolator cores to the
//    target.
//
// Horizon finder array element
//  - Elements near & around each horizon communicate volume data to the same
//    core, sending large buffers over the network.
//  - Volume data is collected on a single core. Possibility of running out of
//    memory on that node.
//  - Interpolation to target points is done in serial on the horizon finder
//    core.
//  - Horizon finder already knows what points to interpolate to. No extra
//    communication needed.
//  - Once sufficient volume data is received, horizon finder can do all fast
//    flow iterations in serial without extra communication.
//
// Overall summary of two methods:
//  - Interpolator method has more charm communication, but does more things in
//    parallel and balances out the load (computation/memory).
//  - Array element method has less charm communication, but puts a larger load
//    (computation/memory) burden on a single core/node
template <size_t Dim>
struct ReceiveInterpolateIterate {
  template <typename ParallelComponent, typename DbTagsList,
            typename Metavariables>
  static void apply(
      db::DataBox<DbTagsList>& box, Parallel::GlobalCache<Metavariables>& cache,
      const int array_index, const Destination destination,
      const size_t global_horizon_find_number, const std::string& frame,
      const LinkedMessageId<double>& temporal_id,
      const ElementId<Dim>& element_id, const Mesh<Dim>& mesh,
      Variables<ah::volume_vars_for_horizon_finding<Dim, Frame::NoFrame>>
          volume_vars_from_element,
      std::vector<std::unique_ptr<Parallel::Callback>> callbacks,
      const bool error_on_failure) {
    const ::Verbosity verbosity = Parallel::get<Tags::Verbosity>(cache);
    const bool quiet_print = verbosity >= ::Verbosity::Quiet;
    const bool verbose_print = verbosity >= ::Verbosity::Verbose;
    const bool debug_print = verbosity >= ::Verbosity::Debug;
    const domain::ObjectLabel object =
        static_cast<domain::ObjectLabel>(array_index);
    const std::string name = "Ah" + domain::name(object);

    const std::set<Storage::NumberAndId>& completed_ids =
        db::get<Tags::CompletedTemporalIds>(box);

    Storage::NumberAndId number_and_id{global_horizon_find_number, temporal_id};

    // This temporal id has finished so we don't need this volume data
    if (completed_ids.count(number_and_id) == 1) {
      return;
    }

    // ReceiveVolumeData
    db::mutate<Tags::VolumeDataAndCallbacks>(
        [&number_and_id, &destination, &frame, &element_id, &mesh,
         &volume_vars_from_element, &callbacks, &error_on_failure](
            const gsl::not_null<std::map<Storage::NumberAndId,
                                         Storage::VolumeDataAndCallbacks>*>
                volume_storage) {
          Storage::VolumeVars volume_vars{mesh,
                                          std::move(volume_vars_from_element)};
          if (UNLIKELY(volume_storage->count(number_and_id) == 0)) {
            volume_storage[number_and_id] = Storage::VolumeDataAndCallbacks{
                frame,
                destination,
                error_on_failure,
                {std::make_pair(element_id, std::move(volume_vars))},
                std::move(callbacks)};
          } else {
            Storage::VolumeDataAndCallbacks& current_vol_data =
                volume_storage->at(number_and_id);
            ASSERT(current_vol_data.frame == frame,
                   "Frames from received volume data don't match.");
            ASSERT(current_vol_data.destination == destination,
                   "Destinations of receive volume data don't match")
            ASSERT(current_vol_data.error_on_failure == error_on_failure,
                   "Received volume data doesn't agree on whether it should "
                   "error on failure.");

            current_vol_data.volume_vars_per_element.emplace(
                element_id, std::move(volume_vars));

            // QUESTION: How do we check these?
            // Only use callbacks from first element? Or make callbacks have
            // some sort of ID?
            current_vol_data.callbacks = std::move(callbacks);
          }
        },
        make_not_null(&box));

    // Check that we have target points for this time and that none of the
    // points are outside the domain. Scoped so we can reuse variable names
    const std::unordered_map<Storage::NumberAndId, Storage::SurfaceAndPoints>&
        all_surface_and_points = db::get<Tags::SurfaceAndPoints>(box);

    // We have not finished a previous horizon find, so we don't yet have
    // points to interpolate to.
    // TODO: What happens if we finish a horizon find, but we don't have any
    // volume data for the next one, but we delete the surface and points for
    // the finished horizon at the end of this while() so we don't have and new
    // points? We could take the previous Strahlkorper in that case? And
    // transform the coords as necessary?
    // Solution: Keep the old number and id around for a while.
    while (all_surface_and_points.count(number_and_id) != 0) {
      const Storage::VolumeDataAndCallbacks& volume_data_and_callbacks =
          db::get<Tags::VolumeDataAndCallbacks>(box).at(number_and_id);
      const Storage::SurfaceAndPoints& surface_and_points =
          all_surface_and_points.at(number_and_id);

      // Do as many iterations as we need
      while (true) {
        // Have to check that all our points are inside the domain
        if (not alg::for_each(surface_and_points.block_logical_coords,
                              [](const auto& optional_block_coords) {
                                return optional_block_coords.has_value();
                              })) {
          ERROR("Apparent horizon finder "
                << name << " failed, reason = "
                << FastFlow::Status::InterpolationFailure
                << ". One or more points is outside of the domain.");
        }

        // Interpolate volume data at this temporal id to the target points that
        // are contained in the elements we already have received
        const bool interpolation_to_all_target_points_complete =
            db::mutate<Tags::InterpolatedVars>(
                [&number_and_id, &volume_data_and_callbacks,
                 &surface_and_points](
                    const gsl::not_null<std::map<Storage::NumberAndId,
                                                 Storage::InterpolatedVars>*>
                        interpolated_vars) {
                  return interpolate_to_target_points(
                      make_not_null(&interpolated_vars->at(number_and_id)),
                      volume_data_and_callbacks.volume_vars_per_element,
                      surface_and_points.block_logical_coords);
                },
                make_not_null(&box));

        // We haven't interpolated to all the target points yet for this
        // iteration. We need to wait for more volume data.
        if (not interpolation_to_all_target_points_complete) {
          return;
        }

        // This edits the Strahlkorper and the previous strahlkorpers in the
        // Storage::SurfaceAndPoints of Tags::SurfaceAndPoints. Also edits the
        // FastFlow
        std::pair<FastFlow::Status, FastFlow::IterInfo> status_and_info =
            single_horizon_find_iteration(make_not_null(&box), number_and_id);

        const FastFlow::Status status = status_and_info.first;
        const FastFlow::IterInfo& info = status_and_info.second;
        const bool horizon_find_has_converged = converged(status);

        if (verbose_print or (quiet_print and horizon_find_has_converged)) {
          Parallel::printf(
              "%s: t=%.16g: its=%d: %.1e<R<%.0e, |R|=%.1g, "
              "|R_grid|=%.1g, %.4g<r<%.4g\n",
              name, number_and_id.id.id, info.iteration, info.min_residual,
              info.max_residual, info.residual_ylm, info.residual_mesh,
              info.r_min, info.r_max);
        }

        // We did a successful iteration. Now we compute the new target points
        // and clear the previous interpolated vars
        if (status == FastFlow::Status::SuccessfulIteration) {
          db::mutate<Tags::SurfaceAndPoints, Tags::InterpolatedVars>(
              [&number_and_id, &cache](
                  const gsl::not_null<std::unordered_map<
                      Storage::NumberAndId, Storage::SurfaceAndPoints>*>
                      local_surface_and_points,
                  const gsl::not_null<std::map<Storage::NumberAndId,
                                               Storage::InterpolatedVars>*>
                      interpolated_vars,
                  const ::FastFlow& fast_flow) {
                Storage::SurfaceAndPoints& surface_and_points =
                    local_surface_and_points->at(number_and_id);

                // Compute the new target points from the old horizon and the
                // fast flow information. Both in NoFrame and the block logical
                // frame
                surface_and_points.cartesian_coords = compute_target_points(
                    surface_and_points.strahlkorper, fast_flow);
                surface_and_points.block_logical_coords =
                    compute_block_logical_target_points(
                        cache, number_and_id.id.id,
                        surface_and_points.cartesian_coords,
                        surface_and_points.frame);

                // Reset the interpolated vars for the next iteration
                interpolated_vars->at(number_and_id) =
                    Storage::InterpolatedVars{};
              },
              make_not_null(&box), db::get<Tags::FastFlow>(box));

          // We need to do another iteration so we go back to the top of this
          // while loop
          continue;
        } else if (not horizon_find_has_converged) {
          // An error occurred during this iteration. Error if we should.
          // Otherwise, clear the interpolated vars and continue on.
          if (volume_data_and_callbacks.error_on_failure) {
            ERROR("Apparent horizon finder " << name
                                             << " failed, reason = " << status);
          } else {
            // Clean everything up for the next horizon find
            clean_up_and_go_to_next_time(make_not_null(&box),
                                         make_not_null(&number_and_id));

            // Break out of the inner while loop and go to the next time if it
            // exists
            break;
          }
        } else {
          // The interpolated variables have been interpolated from the volume
          // to the points on a prolonged_strahlkorper, not to the points on the
          // actual strahlkorper.  So here we do a restriction of these
          // quantities onto the actual strahlkorper.
          db::mutate<Tags::InterpolatedVars>(
              [&number_and_id, &surface_and_points](
                  const gsl::not_null<std::map<Storage::NumberAndId,
                                               Storage::InterpolatedVars>*>
                      interpolated_vars,
                  const FastFlow& fast_flow) {
                using vars_tags =
                    ah::volume_vars_for_horizon_finding<3, Frame::NoFrame>;
                Variables<vars_tags>& current_vars =
                    interpolated_vars->at(number_and_id).vars;

                const ylm::Strahlkorper<Frame::NoFrame>& strahlkorper =
                    surface_and_points.strahlkorper;
                const size_t L_mesh = fast_flow.current_l_mesh(strahlkorper);
                const auto prolonged_strahlkorper =
                    ylm::Strahlkorper<Frame>(L_mesh, L_mesh, strahlkorper);
                auto new_vars = ::Variables<vars_tags>(
                    strahlkorper.ylm_spherepack().physical_size());

                tmpl::for_each<vars_tags>(
                    [&strahlkorper, &prolonged_strahlkorper, &current_vars,
                     &new_vars](auto tag_v) {
                      using tag = typename decltype(tag_v)::type;
                      const auto& old_var = get<tag>(current_vars);
                      auto& new_var = get<tag>(new_vars);
                      auto old_iter = old_var.begin();
                      auto new_iter = new_var.begin();
                      for (; old_iter != old_var.end() and
                             new_iter != new_var.end();
                           ++old_iter, ++new_iter) {
                        *new_iter = strahlkorper.ylm_spherepack().spec_to_phys(
                            prolonged_strahlkorper.ylm_spherepack()
                                .prolong_or_restrict(
                                    prolonged_strahlkorper.ylm_spherepack()
                                        .phys_to_spec(*old_iter),
                                    strahlkorper.ylm_spherepack()));
                      }
                    });
                current_vars = std::move(new_vars);
              },
              make_not_null(&box), db::get<Tags::FastFlow>(box));

          // TODO: Maybe...transform coords to a different frame using
          // strahlkorper_coords_in_different_frame

          // TODO: Interpolate extra vars to target points if we have any.

          // If this was a control system horizon find, update the previous
          // strahlkorpers. We only do this for CS horizon finds because we want
          // the control system logic to be reproducible no matter when we
          // choose to observe the horizon. We also do this before the callbacks
          // in case any of the callbacks need the previous strahlkorpers with
          // the current strahlkorper already in it.
          if (volume_data_and_callbacks.destination ==
              Destination::ControlSystem) {
            db::mutate<Tags::PreviousHorizons>(
                [&number_and_id, &surface_and_points](
                    const gsl::not_null<std::deque<
                        std::pair<double, ylm::Strahlkorper<Frame::NoFrame>>>*>
                        previous_strahlkorpers) {
                  const ylm::Strahlkorper<Frame::NoFrame>& strahlkorper =
                      surface_and_points.strahlkorper;

                  // This is the number of previous strahlkorpers that we
                  // keep around.
                  const size_t num_previous_strahlkorpers = 3;

                  // Save a new previous_strahlkorper.
                  previous_strahlkorpers->emplace_front(number_and_id.id.id,
                                                        strahlkorper);

                  // Remove old previous_strahlkorpers that are no longer
                  // relevant.
                  while (previous_strahlkorpers->size() >
                         num_previous_strahlkorpers) {
                    previous_strahlkorpers->pop_back();
                  }
                },
                make_not_null(&box));
          }

          // Call callbacks
          // QUESTION: These probably shouldn't be parallel callbacks. They
          // should be a new type of callback specific for interpolation that
          // runs on the target. Do these need to be mutable to call them? Or
          // should it not be able to edit anything in the box?
          for (const auto& callback : volume_data_and_callbacks.callbacks) {
            callback->invoke();
          }

          // Clean everything up for the next horizon find
          clean_up_and_go_to_next_time(make_not_null(&box),
                                       make_not_null(&number_and_id));

          // Break out of the inner while loop and go to the next time if it
          // exists
          break;
        }
      }  // while(true)
    }    // while(all_surface_and_points.count(number_and_id) != 0)
  }
};

}  // namespace Actions
}  // namespace ah
