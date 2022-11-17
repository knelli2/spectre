// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <pup.h>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "Domain/Creators/Factory1D.hpp"
#include "Domain/Creators/Factory2D.hpp"
#include "Domain/Creators/Factory3D.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/MinimumGridSpacing.hpp"
#include "Domain/Protocols/Metavariables.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "IO/Observer/Actions/RegisterWithObservers.hpp"
#include "IO/Observer/ArrayComponentId.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Protocols/ReductionDataFormatter.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/VolumeActions.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/Options.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Time/Triggers/SlabCompares.hpp"
#include "Time/Triggers/TimeCompares.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
using MinGridSpacingReductionData = Parallel::ReductionData<
    // Time
    Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
    // Minimum grid spacing
    Parallel::ReductionDatum<double, funcl::Min<>>>;

struct MinGridSpacingFormatter
    : tt::ConformsTo<observers::protocols::ReductionDataFormatter> {
  using reduction_data = MinGridSpacingReductionData;
  std::string operator()(const double time, const double min_grid_spacing) {
    return "Time: " + get_output(time) +
           ", Global inertial minimum grid spacing: " +
           get_output(min_grid_spacing);
  }
  void pup(PUP::er& /*p*/) {}
};

namespace Actions {
template <size_t Dim>
struct ExportCoordinates {
  template <typename ParallelComponent, typename DbTagsList,
            typename ArrayIndex>
  static std::pair<observers::TypeOfObservation, observers::ObservationKey>
  register_info(const db::DataBox<DbTagsList>& /*box*/,
                const ArrayIndex& /*array_index*/) {
    return {observers::TypeOfObservation::Volume,
            observers::ObservationKey("ObserveCoords")};
  }

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const double time = get<Tags::Time>(box);

    const auto& mesh = get<domain::Tags::Mesh<Dim>>(box);
    const auto& inv_jacobian =
        db::get<domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                              Frame::Inertial>>(box);
    const auto& inertial_coordinates =
        db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);
    const auto& grid_coordinates =
        db::get<domain::Tags::Coordinates<Dim, Frame::Grid>>(box);
    const auto deriv_inertial_coordinates =
        partial_derivative(inertial_coordinates, mesh, inv_jacobian);
    // Collect volume data
    // Remove tensor types, only storing individual components
    std::vector<TensorComponent> components;
    components.reserve(1 + 2 * Dim + deriv_inertial_coordinates.size());
    for (size_t d = 0; d < Dim; d++) {
      components.emplace_back("InertialCoordinates_" +
                                  inertial_coordinates.component_name(
                                      inertial_coordinates.get_tensor_index(d)),
                              inertial_coordinates.get(d));
    }
    for (size_t d = 0; d < Dim; d++) {
      components.emplace_back(
          "GridCoordinates_" + grid_coordinates.component_name(
                                   grid_coordinates.get_tensor_index(d)),
          grid_coordinates.get(d));
    }

    for (size_t i = 0; i < deriv_inertial_coordinates.size(); ++i) {
      components.emplace_back(
          "DerivInertialCoordinates_" +
              deriv_inertial_coordinates.component_name(
                  deriv_inertial_coordinates.get_tensor_index(i)),
          deriv_inertial_coordinates[i]);
    }

    // Also output the determinant of the inverse jacobian, which measures
    // the expansion and compression of the grid
    const auto& det_inv_jac = db::get<
        domain::Tags::DetInvJacobian<Frame::ElementLogical, Frame::Inertial>>(
        box);
    components.emplace_back(
        db::tag_name<domain::Tags::DetInvJacobian<Frame::ElementLogical,
                                                  Frame::Inertial>>(),
        get(det_inv_jac));
    // Send data to volume observer
    auto& local_observer = *Parallel::local_branch(
        Parallel::get_parallel_component<observers::Observer<Metavariables>>(
            cache));
    Parallel::simple_action<observers::Actions::ContributeVolumeData>(
        local_observer, observers::ObservationId(time, "ObserveCoords"),
        std::string{"/element_data"},
        observers::ArrayComponentId(
            std::add_pointer_t<ParallelComponent>{nullptr},
            Parallel::ArrayIndex<ElementId<Dim>>(element_id)),
        ElementVolumeData{element_id, std::move(components), mesh});
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

struct FindGlobalMinimumGridSpacing {
  template <typename ParallelComponent, typename DbTagsList,
            typename ArrayIndex>
  static std::pair<observers::TypeOfObservation, observers::ObservationKey>
  register_info(const db::DataBox<DbTagsList>& /*box*/,
                const ArrayIndex& /*array_index*/) {
    return {observers::TypeOfObservation::Reduction,
            observers::ObservationKey("min_grid_spacing")};
  }

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim, typename ActionList, typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const double time = get<Tags::Time>(box);
    const double local_min_grid_spacing =
        get<domain::Tags::MinimumGridSpacing<Dim, Frame::Inertial>>(box);
    auto& local_observer = *Parallel::local_branch(
        Parallel::get_parallel_component<observers::Observer<Metavariables>>(
            cache));
    Parallel::simple_action<observers::Actions::ContributeReductionData>(
        local_observer, observers::ObservationId(time, "min_grid_spacing"),
        observers::ArrayComponentId{
            std::add_pointer_t<ParallelComponent>{nullptr},
            Parallel::ArrayIndex<ElementId<Dim>>(element_id)},
        std::string{"/MinGridSpacing"},
        std::vector<std::string>{"Time", "MinGridSpacing"},
        MinGridSpacingReductionData{time, local_min_grid_spacing},
        std::make_optional(MinGridSpacingFormatter{}));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace Actions

template <size_t Dim, bool EnableTimeDependentMaps>
struct Metavariables {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool local_time_stepping = false;
  // A placeholder system for the domain creators
  struct system {};

  struct domain : tt::ConformsTo<::domain::protocols::Metavariables> {
    static constexpr bool enable_time_dependent_maps = EnableTimeDependentMaps;
  };

  using const_global_cache_tags = tmpl::list<>;

  static constexpr Options::String help{
      "Export the inertial coordinates of the Domain specified in the input "
      "file. The output can be used to compute initial data externally, for "
      "instance. Also outputs the determinant of the inverse jacobian as a "
      "diagnostic of Domain quality: values far from unity indicate "
      "compression or expansion of the grid."};

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<StepChooser<StepChooserUse::LtsStep>, tmpl::list<>>,
        tmpl::pair<StepChooser<StepChooserUse::Slab>, tmpl::list<>>,
        tmpl::pair<Event, tmpl::list<Events::Completion>>,
        tmpl::pair<TimeStepper, TimeSteppers::time_steppers>,
        tmpl::pair<Trigger,
                   tmpl::list<Triggers::SlabCompares, Triggers::TimeCompares>>>;
  };

  using component_list = tmpl::list<
      DgElementArray<
          Metavariables,
          tmpl::list<
              Parallel::PhaseActions<
                  Parallel::Phase::Initialization,
                  tmpl::list<
                      Initialization::Actions::TimeAndTimeStep<Metavariables>,
                      evolution::dg::Initialization::Domain<Dim>,
                      Initialization::Actions::AddComputeTags<
                          ::domain::Tags::MinimumGridSpacingCompute<
                              Dim, Frame::Inertial>>,
                      ::Initialization::Actions::
                          RemoveOptionsAndTerminatePhase>>,
              Parallel::PhaseActions<
                  Parallel::Phase::Register,
                  tmpl::list<observers::Actions::RegisterWithObservers<
                                 Actions::ExportCoordinates<Dim>>,
                             observers::Actions::RegisterWithObservers<
                                 Actions::FindGlobalMinimumGridSpacing>,
                             Parallel::Actions::TerminatePhase>>,
              Parallel::PhaseActions<
                  Parallel::Phase::Execute,
                  tmpl::list<Actions::AdvanceTime,
                             Actions::ExportCoordinates<Dim>,
                             Actions::FindGlobalMinimumGridSpacing,
                             Actions::RunEventsAndTriggers>>>>,
      observers::Observer<Metavariables>,
      observers::ObserverWriter<Metavariables>>;

  using observed_reduction_data_tags = observers::make_reduction_data_tags<
      tmpl::list<MinGridSpacingReductionData>>;

  static constexpr std::array<Parallel::Phase, 4> default_phase_order{
      {Parallel::Phase::Initialization, Parallel::Phase::Register,
       Parallel::Phase::Execute, Parallel::Phase::Exit}};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling,
    &setup_memory_allocation_failure_reporting,
    &disable_openblas_multithreading,
    &domain::creators::register_derived_with_charm,
    &domain::creators::time_dependence::register_derived_with_charm,
    &domain::FunctionsOfTime::register_derived_with_charm,
    &Parallel::register_factory_classes_with_charm<metavariables>};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
/// \endcond
