// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <ostream>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/Systems/Cce/Actions/Psi0Matching.hpp"
#include "Evolution/Systems/Cce/ReceiveTags.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdHelpers.hpp"

namespace logging::Tags {
template <typename OptionsGroup>
struct Verbosity;
}

namespace Cce::Actions {
struct ReceiveGhNextTime {
  // CCE only works in 3D so don't bother templating this on the volume
  // dimension
  using inbox_tags = tmpl::list<Cce::ReceiveTags::GhNextTimeToCcm>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*component*/) {
    auto& inbox = tuples::get<Cce::ReceiveTags::GhNextTimeToCcm>(inboxes);

    const ::Verbosity verbosity =
        Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(cache);
    const size_t expected_number_of_boundary_elements =
        db::get<Tags::ExpectedNumberOfBoundaryElements>(box);
    const TimeStepId& current_ccm_time = db::get<::Tags::TimeStepId>(box);
    const TimeStepId& next_ccm_time =
        db::get<::Tags::Next<::Tags::TimeStepId>>(box);

    // We have not received a next time from all the elements, so wait until we
    // have. This allows us to potentially do less dense output if we have
    // duplicate next GH times.
    if (inbox.size() != expected_number_of_boundary_elements) {
      // Debug because this will print for every GH boundary element
      if (verbosity >= ::Verbosity::Debug) {
        Parallel::printf(
            "ReceiveGhNextTime, t = %.16f, next_t = %.16f: Received only %d "
            "next GH times from %d total boundary elements. Waiting for the "
            "rest.\n",
            current_ccm_time.substep_time(), next_ccm_time.substep_time(),
            inbox.size(), expected_number_of_boundary_elements);
      }
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    } else {
      if (verbosity >= ::Verbosity::Verbose) {
        Parallel::printf(
            "ReceiveGhNextTime, t = %.16f, next_t = %.16f: Received next times "
            "from all %d GH boundary elements\n",
            current_ccm_time.substep_time(), next_ccm_time.substep_time(),
            expected_number_of_boundary_elements);
      }
    }

    // We now have a time for every element. (Potentially) Do some dense output
    std::map<TimeStepId, std::vector<ElementId<3>>> times_for_dense_output{};
    // We cannot erase elements as we are looping through the inbox because the
    // iterators to erased elements become invalidated and we get UB.
    std::unordered_set<ElementId<3>> elements_to_erase{};

    // Loop over all elements and gather the ones we need to do dene output to.
    for (const auto& [element, times] : inbox) {
      const TimeStepId& next_gh_time = times.second;

      if (next_gh_time.substep_time() >= next_ccm_time.substep_time()) {
        // If the next time for this element is after our next time, we don't
        // need to do dense output to this element so we can skip
        continue;
      } else {
        // The next time for this element is before our next time so we have
        // to do dense output. We aggregate the times we need to do dense
        // output to so we don't do more work than necessary.
        if (times_for_dense_output.find(next_gh_time) !=
            times_for_dense_output.end()) {
          times_for_dense_output.at(next_gh_time).emplace_back(element);
        } else {
          times_for_dense_output[next_gh_time] =
              std::vector<ElementId<3>>{element};
        }

        // Since this element needed to do dense output, we no longer keep it in
        // the inbox. We'll have to wait again to see if the next time it sends
        // is after our next time.
        elements_to_erase.insert(element);
      }
    }

    for (const ElementId<3>& element : elements_to_erase) {
      inbox.erase(element);
    }

    // The only way we are allowed to continue is if no elements need to do
    // dense output. This signifies that all elements have a next time that is
    // after our next time so we must continue on.
    if (times_for_dense_output.empty()) {
      if (verbosity >= ::Verbosity::Verbose) {
        Parallel::printf(
            "ReceiveGhNextTime, t = %.16f, next_t = %.16f: All received next "
            "GH times are after current CCE next time. Continuing with CCE "
            "evolution.\n",
            current_ccm_time.substep_time(), next_ccm_time.substep_time());
      }

      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    } else if (verbosity >= ::Verbosity::Verbose) {
      Parallel::printf(
          "ReceiveGhNextTime, t = %.16f, next_t = %.16f: Need to do dense "
          "output to the following times: %s\n",
          current_ccm_time.substep_time(), next_ccm_time.substep_time(),
          keys_of(
              times_for_dense_output,
              [](std::ostream& out,
                 typename std::map<TimeStepId, std::vector<ElementId<3>>>::
                     const_iterator it) { out << it->first.substep_time(); }));
    }

    // Calculate Psi0!
    // For each dense_output_time_id, we will need to
    //  1. Reset BondiJ to its value at current_ccm_time
    //  1. Get BondiJ at dense_output_time_id by doing dense output
    //  2. Calculate Psi0
    //  3. Send Psi0 to the elements that need it at that time
    using bondi_j_var_tag = ::Tags::Variables<tmpl::list<Tags::BondiJ>>;
    using temp_bondi_j_var_tag =
        ::Tags::Variables<tmpl::list<Tags::TemporaryBondiJ>>;
    using history_tag = ::Tags::HistoryEvolvedVariables<bondi_j_var_tag>;
    using psi0_mutators =
        typename CalculatePsi0AndDerivAtInnerBoundary::mutators;

    db::mutate<temp_bondi_j_var_tag>(
        make_not_null(&box),
        [](const gsl::not_null<typename temp_bondi_j_var_tag::type*>
               temp_bondi_j,
           const typename bondi_j_var_tag::type& bondi_j) {
          temp_bondi_j->initialize(bondi_j.number_of_grid_points());
        },
        db::get<bondi_j_var_tag>(box));

    const auto set_variables = [&box](const auto var_tag_to_set_v,
                                      const auto var_tag_to_set_from_v) {
      using var_tag_to_set = std::decay_t<decltype(var_tag_to_set_v)>;
      using var_tag_to_set_from = std::decay_t<decltype(var_tag_to_set_from_v)>;
      using tag_to_set = tmpl::front<typename var_tag_to_set::tags_list>;
      using tag_to_set_from =
          tmpl::front<typename var_tag_to_set_from::tags_list>;
      db::mutate<var_tag_to_set>(
          make_not_null(&box),
          [](const gsl::not_null<typename var_tag_to_set::type*>
                 temporary_bondi_j,
             const typename var_tag_to_set_from::type& current_bondi_j) {
            get<tag_to_set>(*temporary_bondi_j) =
                get<tag_to_set_from>(current_bondi_j);
          },
          db::get<var_tag_to_set_from>(box));
    };

    // Save the current value of BondiJ so we can reset it afterwards
    set_variables(temp_bondi_j_var_tag{}, bondi_j_var_tag{});

    // Loop over all times to do dense output to and send the result to all
    // elements waiting at that time
    bool first_dense_output = true;
    for (const auto& time_and_elements : times_for_dense_output) {
      const TimeStepId& dense_output_time_id = time_and_elements.first;
      const std::vector<ElementId<3>>& elements = time_and_elements.second;
      // Reset BondiJ. If this is the first dense output time, save ourselves a
      // copy because BondiJ is already at the correct values
      if (not first_dense_output) {
        set_variables(bondi_j_var_tag{}, temp_bondi_j_var_tag{});
      } else {
        first_dense_output = false;
      }

      // Actually do the dense output for BondiJ to the dense output time
      db::mutate<bondi_j_var_tag>(
          make_not_null(&box),
          [&dense_output_time_id](
              const gsl::not_null<typename bondi_j_var_tag::type*> bondi_j,
              const typename history_tag::type& history,
              const auto& time_stepper) {
            time_stepper.dense_update_u(bondi_j, history,
                                        dense_output_time_id.substep_time());
          },
          db::get<history_tag>(box), db::get<::Tags::TimeStepper<>>(box));

      // Calculate Psi0 on the inner boundary (stored as
      // ::Tags::Variables<metavars::ccm_psi0>)
      tmpl::for_each<psi0_mutators>([&box](auto mutator_v) {
        using mutator = typename decltype(mutator_v)::type;
        db::mutate_apply<mutator>(make_not_null(&box));
      });

      // Once we have Psi0, send it to all the boundary elements
      for (const ElementId<3>& element : elements) {
        Parallel::receive_data<
            Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
            Parallel::get_parallel_component<
                typename Metavariables::gh_dg_element_array>(cache)[element],
            dense_output_time_id,
            db::get<::Tags::Variables<typename Metavariables::ccm_psi0>>(box),
            false);
      }
    }

    // Restore the BondiJ in the box to what it was before we did dense output
    set_variables(bondi_j_var_tag{}, temp_bondi_j_var_tag{});

    if (verbosity >= ::Verbosity::Verbose) {
      Parallel::printf(
          "ReceiveGhNextTime, t = %.16f, next_t = %.16f: Completed dense "
          "output to all GH elements. Waiting until we are certain all GH "
          "elements' next time is after CCE's next time.\n",
          current_ccm_time.substep_time(), next_ccm_time.substep_time());
    }

    // If we get here, we've sent Psi0 to all the elements we have in our inbox
    // that needed to do dense output. Thus we must wait for new next GH times
    // from all these elements to check whether they are before or after our
    // next time.
    return {Parallel::AlgorithmExecution::Retry, std::nullopt};
  }
};
}  // namespace Cce::Actions
