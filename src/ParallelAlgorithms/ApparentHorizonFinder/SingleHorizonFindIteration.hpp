// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <deque>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace ah {

/// TODO: Update docs
/// \brief post interpolation callback (see
/// intrp::protocols::PostInterpolationCallback) that does a FastFlow iteration
/// and triggers another one until convergence.
///
/// Assumes that InterpolationTargetTag contains an additional
/// type alias called `post_horizon_find_callbacks`, which is a list of
/// structs, each of which has a function
///
/// \snippet
/// ParallelAlgorithms/ApparentHorizonFinder/Test_ApparentHorizonFinder.cpp
/// post_horizon_find_callback_example
///
/// that is called if the FastFlow iteration has converged.
/// InterpolationTargetTag also is assumed to contain an additional
/// struct called `horizon_find_failure_callback`, which has a function
///
/// \snippet
/// ParallelAlgorithms/ApparentHorizonFinder/Test_ApparentHorizonFinder.cpp
/// horizon_find_failure_callback_example
///
/// that is called if the FastFlow iteration or the interpolation has
/// failed.
///
/// Uses:
/// - Metavariables:
///   - `temporal_id`
/// - DataBox:
///   - `logging::Tags::Verbosity<InterpolationTargetTag>`
///   - `::gr::Tags::InverseSpatialMetric<DataVector, 3, Frame>`
///   - `::gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame>`
///   - `::gr::Tags::SpatialChristoffelSecondKind<DataVector, 3, Frame>`
///   - `::ah::Tags::FastFlow`
///   - `ylm::Tags::Strahlkorper<Frame>`
///
/// Modifies:
/// - DataBox:
///   - `::ah::Tags::FastFlow`
///   - `ylm::Tags::Strahlkorper<Frame>`
///
/// This is an InterpolationTargetTag::post_interpolation_callback;
/// see intrp::protocols::InterpolationTargetTag for details on
/// InterpolationTargetTag.
///
/// ### Output
///
/// Optionally, a single line of output is printed to stdout either on each
/// iteration (if verbosity > Verbosity::Quiet) or on convergence
/// (if verbosity > Verbosity::Silent).  The output consists of the
/// following labels, and values associated with each label:
///  - t        = time given by value of `temporal_id` argument to `apply`
///  - its      = current iteration number.
///  - R        = min and max of residual over all prolonged grid points.
///  - |R|      = L2 norm of residual, counting only L modes solved for.
///  - |R_mesh| = L2 norm of residual over prolonged grid points.
///  - r        = min and max radius of trial horizon surface.
///
/// #### Difference between |R| and |R_mesh|:
///  The horizon is represented in a \f$Y_{lm}\f$ expansion up
///  to \f$l=l_{\mathrm{surface}}\f$;
///  the residual |R| represents the failure of that surface to satisfy
///  the apparent horizon equation.
///
///  However, at each iteration we also interpolate the horizon surface
///  to a higher resolution ("prolongation").  The prolonged surface
///  includes \f$Y_{lm}\f$ coefficients up to \f$l=l_{\mathrm{mesh}}\f$,
///  where \f$l_{\mathrm{mesh}} > l_{\mathrm{surface}}\f$.
///  The residual computed on this higher-resolution surface is |R_mesh|.
///
///  As iterations proceed, |R| should decrease until it reaches
///  numerical roundoff error, because we are varying all the \f$Y_{lm}\f$
///  coefficients up to  \f$l=l_{\mathrm{surface}}\f$
///  to minimize the residual.  However,
///  |R_mesh| should eventually stop decreasing as iterations proceed,
///  hitting a floor that represents the numerical truncation error of
///  the solution.
///
///  The convergence criterion looks at both |R| and |R_mesh|:  Once
///  |R| is small enough and |R_mesh| has stabilized, it is pointless
///  to keep iterating (even though one could iterate until |R| reaches
///  roundoff).
///
///  Problems with convergence of the apparent horizon finder can often
///  be diagnosed by looking at the behavior of |R| and |R_mesh| as a
///  function of iteration.
///
template <typename DbTags>
std::pair<FastFlow::Status, FastFlow::IterInfo> single_horizon_find_iteration(
    const gsl::not_null<db::DataBox<DbTags>*> box,
    const Storage::NumberAndId& number_and_id) {
  if (get<Tags::FastFlow>(*box).current_iteration() == 0) {
    // If we get here, we are in a new apparent horizon search, as
    // opposed to a subsequent iteration of the same horizon search.
    //
    // So put new initial guess into ylm::Tags::Strahlkorper<Frame>.
    // We need to do this now, and not at the end of the previous horizon
    // search, because only now do we know the temporal_id of this horizon
    // search.
    db::mutate<Tags::SurfaceAndPoints>(
        [&number_and_id](
            const gsl::not_null<std::unordered_map<Storage::NumberAndId,
                                                   Storage::SurfaceAndPoints>*>
                surface_and_points,
            const std::deque<
                std::pair<double, ylm::Strahlkorper<Frame::NoFrame>>>
                previous_strahlkorpers) {
          Storage::SurfaceAndPoints& current_surface_and_points =
              surface_and_points->at(number_and_id);
          ylm::Strahlkorper<Frame::NoFrame>& strahlkorper =
              current_surface_and_points.strahlkorper;
          // If we have zero previous_strahlkorpers, then the
          // initial guess is already in strahlkorper, so do
          // nothing.
          //
          // If we have one previous_strahlkorper, then we have had
          // a successful horizon find, and the initial guess for the
          // next horizon find is already in strahlkorper, so
          // again we do nothing.
          //
          // If we have 2 previous_strahlkorpers and the time of the second
          // one is a NaN, this means that the corresponding
          // previous_strahlkorper is the original initial guess, so
          // again we do nothing.
          //
          // If we have 2 or more valid previous_strahlkorpers, then
          // we set the initial guess by linear extrapolation in time
          // using the last 2 previous_strahlkorpers.
          if (previous_strahlkorpers->size() > 1 and
              not std::isnan((*previous_strahlkorpers)[1].first)) {
            const double dt_0 =
                (*previous_strahlkorpers)[0].first - number_and_id.id.id;
            const double dt_1 =
                (*previous_strahlkorpers)[1].first - number_and_id.id.id;
            const double fac_0 = dt_1 / (dt_1 - dt_0);
            const double fac_1 = 1.0 - fac_0;
            // Here we assume that
            // * Expansion center of all the Strahlkorpers are equal.
            // * Maximum L of all the Strahlkorpers are equal.
            // It is easy to relax the max L assumption once we start
            // adaptively changing the L of the strahlkorpers.
            strahlkorper->coefficients() =
                fac_0 * (*previous_strahlkorpers)[0].second.coefficients() +
                fac_1 * (*previous_strahlkorpers)[1].second.coefficients();
          }
        },
        box, db::get<Tags::PreviousHorizons>(*box));
  }

  // Get the interpolated vars on the target set of points
  const Variables<ah::volume_vars_for_horizon_finding<3, Frame::NoFrame>>&
      interpolated_vars =
          db::get<Tags::InterpolatedVars>(*box).at(number_and_id).vars;
  const auto& inv_g =
      db::get<::gr::Tags::InverseSpatialMetric<DataVector, 3, Frame::NoFrame>>(
          interpolated_vars);
  const auto& ex_curv =
      db::get<::gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame::NoFrame>>(
          interpolated_vars);
  const auto& christoffel = db::get<
      ::gr::Tags::SpatialChristoffelSecondKind<DataVector, 3, Frame::NoFrame>>(
      interpolated_vars);

  std::pair<FastFlow::Status, FastFlow::IterInfo> status_and_info;

  // Do a FastFlow iteration.
  db::mutate<Tags::FastFlow, Tags::SurfaceAndPoints>(
      [&number_and_id, &inv_g, &ex_curv, &christoffel, &status_and_info](
          const gsl::not_null<::FastFlow*> fast_flow,
          const gsl::not_null<std::unordered_map<Storage::NumberAndId,
                                                 Storage::SurfaceAndPoints>*>
              surface_and_points) {
        auto& strahlkorper = surface_and_points->at(number_and_id).strahlkorper;
        status_and_info = fast_flow->template iterate_horizon_finder<Frame>(
            make_not_null(&strahlkorper), inv_g, ex_curv, christoffel);
      },
      box);

  return status_and_info;
}
}  // namespace ah
