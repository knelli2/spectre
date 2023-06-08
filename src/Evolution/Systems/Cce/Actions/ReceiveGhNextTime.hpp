// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
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
template <bool DoDenseOutput, typename DbTags, typename Metavariables>
void calculate_and_send_psi0(
    const gsl::not_null<db::DataBox<DbTags>*> box,
    Parallel::GlobalCache<Metavariables>& cache,
    const std::map<TimeStepId, std::optional<std::vector<ElementId<3>>>>&
        times_for_output) {
  if (times_for_output.empty()) {
    return;
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
  using psi0_mutators = typename CalculatePsi0AndDerivAtInnerBoundary::mutators;

  const auto set_variables = [&box](const auto var_tag_to_set_v,
                                    const auto var_tag_to_set_from_v) {
    using var_tag_to_set = std::decay_t<decltype(var_tag_to_set_v)>;
    using var_tag_to_set_from = std::decay_t<decltype(var_tag_to_set_from_v)>;
    using tag_to_set = tmpl::front<typename var_tag_to_set::tags_list>;
    using tag_to_set_from =
        tmpl::front<typename var_tag_to_set_from::tags_list>;
    db::mutate<var_tag_to_set>(
        box,
        [](const gsl::not_null<typename var_tag_to_set::type*>
               temporary_bondi_j,
           const typename var_tag_to_set_from::type& current_bondi_j) {
          get<tag_to_set>(*temporary_bondi_j) =
              get<tag_to_set_from>(current_bondi_j);
        },
        db::get<var_tag_to_set_from>(*box));
  };

  if constexpr (DoDenseOutput) {
    db::mutate<temp_bondi_j_var_tag>(
        box,
        [](const gsl::not_null<typename temp_bondi_j_var_tag::type*>
               temp_bondi_j,
           const typename bondi_j_var_tag::type& bondi_j) {
          temp_bondi_j->initialize(bondi_j.number_of_grid_points());
        },
        db::get<bondi_j_var_tag>(*box));

    // Save the current value of BondiJ so we can reset it afterwards
    set_variables(temp_bondi_j_var_tag{}, bondi_j_var_tag{});
  }

  // Loop over all times to do dense output to and send the result to all
  // elements waiting at that time
  bool first_dense_output = true;
  for (const auto& time_and_elements : times_for_output) {
    const TimeStepId& dense_output_time_id = time_and_elements.first;
    const std::optional<std::vector<ElementId<3>>>& elements =
        time_and_elements.second;
    if constexpr (DoDenseOutput) {
      // Reset BondiJ. If this is the first dense output time, save ourselves a
      // copy because BondiJ is already at the correct values
      if (not first_dense_output) {
        set_variables(bondi_j_var_tag{}, temp_bondi_j_var_tag{});
      } else {
        first_dense_output = false;
      }

      // Actually do the dense output for BondiJ to the dense output time
      db::mutate<bondi_j_var_tag>(
          box,
          [&dense_output_time_id](
              const gsl::not_null<typename bondi_j_var_tag::type*> bondi_j,
              const typename history_tag::type& history,
              const auto& time_stepper) {
            time_stepper.dense_update_u(bondi_j, history,
                                        dense_output_time_id.substep_time());
          },
          db::get<history_tag>(*box), db::get<::Tags::TimeStepper<>>(*box));
    }

    // Calculate Psi0 on the inner boundary (stored as
    // ::Tags::Variables<metavars::ccm_psi0>)
    tmpl::for_each<psi0_mutators>([&box](auto mutator_v) {
      using mutator = typename decltype(mutator_v)::type;
      db::mutate_apply<mutator>(box);
    });

    // Once we have Psi0, send it to all the boundary elements. If there were
    // elements specified, send to only those. If not, send to all elements
    if (elements.has_value()) {
      for (const ElementId<3>& element : elements.value()) {
        Parallel::receive_data<
            Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
            Parallel::get_parallel_component<
                typename Metavariables::gh_dg_element_array>(cache)[element],
            dense_output_time_id,
            db::get<::Tags::Variables<typename Metavariables::ccm_psi0>>(*box),
            false);
      }
    } else {
      Parallel::receive_data<
          Cce::ReceiveTags::BoundaryData<typename Metavariables::ccm_psi0>>(
          Parallel::get_parallel_component<
              typename Metavariables::gh_dg_element_array>(cache),
          dense_output_time_id,
          db::get<::Tags::Variables<typename Metavariables::ccm_psi0>>(*box),
          false);
    }
  }

  if constexpr (DoDenseOutput) {
    // Restore the BondiJ in the box to what it was before we did dense output
    set_variables(bondi_j_var_tag{}, temp_bondi_j_var_tag{});
  }
}

template <typename Inbox, typename DbTags, typename Metavariables>
bool received_all_gh_times(const Inbox& inbox, const db::DataBox<DbTags>& box,
                           const Parallel::GlobalCache<Metavariables>& cache,
                           const std::string& action_name) {
  const ::Verbosity verbosity =
      Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(cache);
  const size_t expected_number_of_boundary_elements =
      db::get<Tags::ExpectedNumberOfBoundaryElements>(box);
  const TimeStepId& current_ccm_time = db::get<::Tags::TimeStepId>(box);
  const TimeStepId& next_ccm_time =
      db::get<::Tags::Next<::Tags::TimeStepId>>(box);

  // We have not received a next time from all the elements, so wait until
  // we have.
  if (inbox.size() != expected_number_of_boundary_elements) {
    // Debug because this will print for every GH boundary element
    if (verbosity >= ::Verbosity::Debug) {
      Parallel::printf(
          "%s, t = %.16f, next_t = %.16f: Received only %d "
          "next GH times from %d total boundary elements. Waiting for the "
          "rest.\n",
          action_name, current_ccm_time.substep_time(),
          next_ccm_time.substep_time(), inbox.size(),
          expected_number_of_boundary_elements);
    }

    return false;
  } else {
    if (verbosity >= ::Verbosity::Verbose) {
      Parallel::printf(
          "%s, t = %.16f, next_t = %.16f: Received next "
          "times "
          "from all %d GH boundary elements\n",
          action_name, current_ccm_time.substep_time(),
          next_ccm_time.substep_time(), expected_number_of_boundary_elements);
    }

    return true;
  }
}

enum class CriterionForUse { AtCurrentTime, DenseOutputThisStep };

template <typename Inbox>
void get_times_for_output(
    const gsl::not_null<
        std::map<TimeStepId, std::optional<std::vector<ElementId<3>>>>*>
        times_for_output,
    const gsl::not_null<Inbox*> inbox, const TimeStepId& current_ccm_time,
    const TimeStepId& next_ccm_time, const CriterionForUse criterion_for_use) {
  // We cannot erase elements as we are looping through the inbox because
  // the iterators to erased elements become invalidated and we get UB.
  std::unordered_set<ElementId<3>> elements_to_erase{};

  // Loop over all elements and gather the ones we need to do dene output
  // to.
  for (const auto& [element, times] : *inbox) {
    const TimeStepId& next_gh_time = times.second;

    bool criterion_is_met = false;
    if (criterion_for_use == CriterionForUse::AtCurrentTime) {
      criterion_is_met =
          next_gh_time.substep_time() == current_ccm_time.substep_time();
    } else if (criterion_for_use == CriterionForUse::DenseOutputThisStep) {
      criterion_is_met =
          next_gh_time.substep_time() < next_ccm_time.substep_time();
    } else {
      ERROR("Unknown criterion in Cce::Actions::get_times_for_output: "
            << static_cast<int>(criterion_for_use));
    }

    if (criterion_is_met) {
      if (times_for_output->find(next_gh_time) != times_for_output->end()) {
        times_for_output->at(next_gh_time).value().emplace_back(element);
      } else {
        (*times_for_output)[next_gh_time] = std::vector<ElementId<3>>{element};
      }

      // Since we are sending Psi0 to this element now, we no longer keep it
      // in the inbox. We'll have to wait again to see if the next time it
      // sends matches the criterion.
      elements_to_erase.insert(element);
    } else {
      // If the next time for this element doesn't match the criterion, we don't
      // need to remove it so we can skip
      continue;
    }
  }

  for (const ElementId<3>& element : elements_to_erase) {
    inbox->erase(element);
  }
}

template <typename ListOfActionsToSolveForBondiJ>
struct SendPsi0Early {
  using inbox_tags = tmpl::list<Cce::ReceiveTags::GhNextTimeToCcm>;
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList action_list,
      const ParallelComponent* const component) {
    auto& inbox = tuples::get<Cce::ReceiveTags::GhNextTimeToCcm>(inboxes);

    const ::Verbosity verbosity =
        Parallel::get<logging::Tags::Verbosity<Cce::OptionTags::Cce>>(cache);
    const TimeStepId& current_ccm_time = db::get<::Tags::TimeStepId>(box);
    const TimeStepId& next_ccm_time =
        db::get<::Tags::Next<::Tags::TimeStepId>>(box);

    std::map<TimeStepId, std::optional<std::vector<ElementId<3>>>>
        times_for_send{};

    const bool at_slab_zero = current_ccm_time.slab_number() == 0 and
                              current_ccm_time.is_at_slab_boundary();
    if (UNLIKELY(at_slab_zero)) {
      if (verbosity >= ::Verbosity::Verbose) {
        Parallel::printf(
            "SendPsi0Early, t = %.16f, next_t = %.16f: Sending Psi0 early "
            "because we are at slab 0.\n",
            current_ccm_time.substep_time(), next_ccm_time.substep_time());
      }
      tmpl::for_each<ListOfActionsToSolveForBondiJ>(
          [&box, &inboxes, &cache, &array_index, &action_list,
           &component](auto action_v) {
            using action = tmpl::type_from<decltype(action_v)>;
            action::apply(box, inboxes, cache, array_index, action_list,
                          component);
          });

      times_for_send[current_ccm_time] = std::nullopt;

      calculate_and_send_psi0<false>(make_not_null(&box), cache,
                                     times_for_send);

      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    if (not received_all_gh_times(inbox, box, cache, "SendPsi0Early")) {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    get_times_for_output(make_not_null(&times_for_send), make_not_null(&inbox),
                         current_ccm_time, next_ccm_time,
                         CriterionForUse::AtCurrentTime);

    // We don't have to Retry here because we only care about times that are at
    // our current time. So for any times we removed from the inbox because they
    // are at our current time, we are guaranteed that the next time will be
    // after the current time, meaning we don't have to wait. If the next time
    // is before our next time, then we'll need to do dense output in another
    // action. If it's at our next time, then it'll hit this action on the next
    // time through the action list. If it's after our next time, then we don't
    // care either way.
    calculate_and_send_psi0<false>(make_not_null(&box), cache, times_for_send);

    if (verbosity >= ::Verbosity::Verbose) {
      Parallel::printf(
          "SendPsi0Early, t = %.16f, next_t = %.16f: Need to send data to the "
          "following times: %s\n",
          current_ccm_time.substep_time(), next_ccm_time.substep_time(),
          keys_of(
              times_for_send,
              [](std::ostream& out,
                 typename std::map<TimeStepId,
                                   std::optional<std::vector<ElementId<3>>>>::
                     const_iterator it) { out << it->first.substep_time(); }));
    }

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct ReceiveGhNextTime {
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
    const TimeStepId& current_ccm_time = db::get<::Tags::TimeStepId>(box);
    const TimeStepId& next_ccm_time =
        db::get<::Tags::Next<::Tags::TimeStepId>>(box);

    // We have not received a next time from all the elements, so wait until
    // we have. This allows us to potentially do less dense output if we have
    // duplicate next GH times.
    if (not received_all_gh_times(inbox, box, cache, "ReceiveGhNextTime")) {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    // We now have a time for every element. (Potentially) Do some dense
    // output
    std::map<TimeStepId, std::optional<std::vector<ElementId<3>>>>
        times_for_dense_output{};

    get_times_for_output(make_not_null(&times_for_dense_output),
                         make_not_null(&inbox), current_ccm_time, next_ccm_time,
                         CriterionForUse::DenseOutputThisStep);

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
                 typename std::map<TimeStepId,
                                   std::optional<std::vector<ElementId<3>>>>::
                     const_iterator it) { out << it->first.substep_time(); }));
    }

    calculate_and_send_psi0<true>(make_not_null(&box), cache,
                                  times_for_dense_output);

    if (verbosity >= ::Verbosity::Verbose) {
      Parallel::printf(
          "ReceiveGhNextTime, t = %.16f, next_t = %.16f: Completed dense "
          "output to all GH elements. Waiting until we are certain all GH "
          "elements' next time is after CCE's next time.\n",
          current_ccm_time.substep_time(), next_ccm_time.substep_time());
    }

    // If we get here, we've sent Psi0 to all the elements we have in our
    // inbox that needed to do dense output. Thus we must wait for new next GH
    // times from all these elements to check whether they are before or after
    // our next time.
    return {Parallel::AlgorithmExecution::Retry, std::nullopt};
  }
};
}  // namespace Cce::Actions
