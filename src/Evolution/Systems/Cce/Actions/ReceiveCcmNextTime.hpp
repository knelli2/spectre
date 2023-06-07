// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/Systems/Cce/IsBoundaryElement.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/Gsl.hpp"

#include "Parallel/Printf.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"

namespace logging::Tags {
template <typename OptionsGroup>
struct Verbosity;
}

namespace Cce {
namespace Actions {

template <typename WorldtubeTargetTag>
struct ReceiveCcmNextTime {
  // CCE only works in 3D so don't bother templating this on the volume
  // dimension
  using inbox_tags =
      tmpl::list<Cce::ReceiveTags::CcmNextTimeToGH,
                 evolution::dg::Tags::BoundaryCorrectionAndGhostCellsInbox<3>>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const component) {
    auto& inbox = tuples::get<Cce::ReceiveTags::CcmNextTimeToGH>(inboxes);

    const double gh_time = db::get<::Tags::TimeStepId>(box).substep_time();

    const bool debug_print =
        Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(cache) >=
        ::Verbosity::Debug;

    if (not Cce::is_outer_boundary_element<true>(make_not_null(&inbox), box,
                                                 cache)) {
      // if (debug_print) {
      //   Parallel::printf(
      //       "ReceiveCcmNextTime, %s: Not a boundary element. Continuing.\n",
      //       gh_time, array_index);
      // }
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    // If we don't have a next time for CCE yet, wait until we have one. We
    // can't continue on because we don't know if the next CCE time will be
    // before the next GH time (meaning we need to do dense output) or after
    // (meaning we can continue on). We do this after we know we are an outer
    // boundary element so we don't hold up volume elements.
    if (inbox.size() == 0) {
      if (debug_print) {
        Parallel::printf(
            "ReceiveCcmNextTime, %s: Current gh time is %.16f. No times in "
            "inbox, waiting\n",
            get_output(array_index), gh_time);
      }
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    // Now loop over the times in the inbox. We do this because there may be
    // multiple next CCE times we have to do dense output to before the next GH
    // time.
    while (inbox.size() > 0) {
      const TimeStepId& current_gh_time = db::get<::Tags::TimeStepId>(box);
      const TimeStepId& next_gh_time =
          db::get<::Tags::Next<::Tags::TimeStepId>>(box);
      const TimeStepId& next_cce_time = inbox.begin()->second;

      // Use CCE next time as scale because that will always be equal to or
      // larger than current GH time
      const bool cce_next_time_at_current_gh_time =
          current_gh_time.substep_time() == next_cce_time.substep_time();
      const bool cce_next_time_after_next_gh_time =
          next_cce_time.substep_time() > next_gh_time.substep_time();

      if (debug_print) {
        Parallel::printf(
            "ReceiveCcmNextTime, %s: Current gh time is %.16f. Inbox has %d "
            "times. Next gh time is %.16f, next CCE time is %.16f. We are%s "
            "doing dense output this step.\n",
            get_output(array_index), gh_time, inbox.size(),
            next_gh_time.substep_time(), next_cce_time.substep_time(),
            (cce_next_time_after_next_gh_time ? " not" : ""));
      }

      // If the first next CCE time is after (or at) the next GH time, continue
      // on because we don't need to do dense output this GH time step. Don't
      // clear the inbox because the next time through the action list we'll
      // need it.
      if (cce_next_time_after_next_gh_time) {
        return {Parallel::AlgorithmExecution::Continue, std::nullopt};
      }

      // We have to have all boundary corrections before we can do dense output.
      // We do this now so that we aren't unnecessarily waiting for boundary
      // corrections if the next CCE time is after our current time step
      if (not cce_next_time_at_current_gh_time) {
        const double current_time = db::get<::Tags::Time>(box);

        const auto set_time_to = [&box](const double new_time) {
          db::mutate<::Tags::Time>(
              make_not_null(&box),
              [&new_time](const gsl::not_null<double*> time) {
                *time = new_time;
              });
        };

        // The receive_boundary_data_lts function uses the value of ::Tags::Time
        // to represent the dense output time, so we need to set it to the
        // correct time before we call the function and then set it back after.
        set_time_to(next_cce_time.substep_time());

        // DENSE OUTPUT!!!!!!!!!!!!!!!!!!!!!!!1
        if (not evolution::dg::receive_boundary_data_local_time_stepping<
                Metavariables, true>(make_not_null(&box),
                                     make_not_null(&inboxes))) {
          // set_time_to(current_time);
          if (debug_print) {
            Parallel::printf(
                "ReceiveCcmNextTime, %s: Current gh time is %.16f. Waiting for "
                "mortar data. Next cce time is %f\n",
                get_output(array_index), gh_time, next_cce_time.substep_time());
          }
          return {Parallel::AlgorithmExecution::Retry, std::nullopt};
        }

        set_time_to(current_time);
      }

      using system = gh::System<3>;
      using variables_tag = typename system::variables_tag;
      using evolved_vars_list = typename variables_tag::tags_list;
      using variables_of_evolved_vars = typename variables_tag::type;

      // The first next CCE time is before the next GH time which means we need
      // to do dense output
      const auto& time_stepper = db::get<::Tags::TimeStepper<>>(box);
      const auto& history_evolved_vars =
          db::get<::Tags::HistoryEvolvedVariables<variables_tag>>(box);
      const auto& variables_evolved_vars = db::get<variables_tag>(box);

      variables_of_evolved_vars evolved_vars = variables_evolved_vars;

      // Even if we are at the correct time, since ComputeTimeDerivative has
      // already been called, we have to do a dense update. Otherwise the
      // variables are not in the state we need them to be
      bool succeeded = time_stepper.dense_update_u(
          make_not_null(&evolved_vars), history_evolved_vars,
          next_cce_time.substep_time());

      if (not cce_next_time_at_current_gh_time) {
        // Apply boundary corrections. We don't use db::apply because we need to
        // pass in our own variables and dense output time.
        evolution::dg::ApplyBoundaryCorrections<true, system, 3, true>::apply(
            make_not_null(&evolved_vars),
            db::get<evolution::dg::Tags::MortarDataHistory<
                3, db::add_tag_prefix<::Tags::dt, variables_tag>::type>>(box),
            db::get<domain::Tags::Mesh<3>>(box),
            db::get<evolution::dg::Tags::MortarMesh<3>>(box),
            db::get<evolution::dg::Tags::MortarSize<3>>(box),
            db::get<::dg::Tags::Formulation>(box),
            db::get<evolution::dg::Tags::NormalCovectorAndMagnitude<3>>(box),
            db::get<::Tags::TimeStepper<>>(box),
            db::get<evolution::Tags::BoundaryCorrection<system>>(box),
            next_cce_time.substep_time());

        if (debug_print) {
          Parallel::printf(
              "ReceiveCcmNextTime, %s: Current gh time is %.16f. Dense output "
              "to %f%s complete\n",
              get_output(array_index), gh_time, next_cce_time.substep_time(),
              (succeeded ? "" : " not"));
          // Parallel::printf("Vars:\n%s\n\n", evolved_vars);
        }
      } else {
        if (debug_print) {
          Parallel::printf(
              "ReceiveCcmNextTime, %s: Current gh time is %.16f. No need to do "
              "dense output to %f because we are already at that time\n",
              get_output(array_index), gh_time, next_cce_time.substep_time());
        }
      }

      auto interpolate_event =
          ::intrp::Events::InterpolateWithoutInterpComponent<
              3, WorldtubeTargetTag, Metavariables, evolved_vars_list>{};

      // These arguments are needed to pass to the interpolate_event
      const auto& point_info =
          db::get<intrp::Tags::InterpPointInfo<Metavariables>>(box);
      const auto& mesh = db::get<domain::Tags::Mesh<3>>(box);

      // These must be in this order when passed to the interpolating event
      // because that's the order they are in the `evolved_vars_list` type alias
      const auto& spacetime_metric =
          get<gr::Tags::SpacetimeMetric<DataVector, 3, ::Frame::Inertial>>(
              evolved_vars);
      const auto& pi =
          get<gh::Tags::Pi<DataVector, 3, ::Frame::Inertial>>(evolved_vars);
      const auto& phi =
          get<gh::Tags::Phi<DataVector, 3, ::Frame::Inertial>>(evolved_vars);

      // Actually do the interpolation to the worldtube target. The TimeStepId
      // associated with current inbox value is the TimeStepId we need to pass
      // to the interpolation target, and then to CCE because CCE is expecting
      // this particular TimeStepId. We call the operator() by hand here rather
      // than use a ::apply() call because we don't want to pass the evolved
      // variables from the box to this interpolation. They are at the incorrect
      // time. We instead want to pass out dense-outputted variables which are
      // at the next CCE time.
      interpolate_event(next_cce_time, point_info, mesh, spacetime_metric, pi,
                        phi, cache, array_index, component);

      // Now that we have sent GH data at the next CCE time, remove it from the
      // inbox
      inbox.erase(inbox.begin());
    }

    if (debug_print) {
      Parallel::printf(
          "ReceiveCcmNextTime, %s: Current gh time is %.16f. Finished with all "
          "CCE times. Waiting\n",
          get_output(array_index), gh_time);
    }

    // If we get here, we've sent GH data to all the next CCE times we have in
    // our inbox and all of those next CCE times happened to be before our next
    // GH time. Thus we must wait for a new next CCE time to check whether it is
    // before or after the next GH time.
    return {Parallel::AlgorithmExecution::Retry, std::nullopt};
  }
};

}  // namespace Actions
}  // namespace Cce
