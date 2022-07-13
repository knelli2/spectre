// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstdint>
#include <vector>

#include "ApparentHorizons/ComputeHorizonVolumeQuantities.hpp"
#include "ApparentHorizons/ComputeHorizonVolumeQuantities.tpp"
#include "ApparentHorizons/ComputeItems.hpp"
#include "ApparentHorizons/HorizonAliases.hpp"
#include "ApparentHorizons/ObserveCenters.hpp"
#include "ApparentHorizons/Tags.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Creators/Factory1D.hpp"
#include "Domain/Creators/Factory2D.hpp"
#include "Domain/Creators/Factory3D.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/FunctionsOfTimeAreReady.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/Protocols/Metavariables.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsCharacteristicSpeeds.hpp"
#include "Evolution/Actions/RunEventsAndDenseTriggers.hpp"
#include "Evolution/ComputeTags.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "Evolution/EventsAndDenseTriggers/DenseTriggers/Factory.hpp"
#include "Evolution/Executables/Cce/CharacteristicExtractBase.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Initialization/NonconservativeSystem.hpp"
#include "Evolution/Initialization/SetVariables.hpp"
#include "Evolution/NumericInitialData.hpp"
#include "Evolution/Systems/Cce/Actions/InterpolateDuringSelfStart.hpp"
#include "Evolution/Systems/Cce/BoundaryData.hpp"
#include "Evolution/Systems/Cce/Components/CharacteristicEvolution.hpp"
#include "Evolution/Systems/Cce/Components/WorldtubeBoundary.hpp"
#include "Evolution/Systems/Cce/Initialize/ConformalFactor.hpp"
#include "Evolution/Systems/Cce/Initialize/InitializeJ.hpp"
#include "Evolution/Systems/Cce/Initialize/InverseCubic.hpp"
#include "Evolution/Systems/Cce/Initialize/NoIncomingRadiation.hpp"
#include "Evolution/Systems/Cce/Initialize/RegisterInitializeJWithCharm.hpp"
#include "Evolution/Systems/Cce/Initialize/ZeroNonSmooth.hpp"
#include "Evolution/Systems/Cce/IntegrandInputSteps.hpp"
#include "Evolution/Systems/Cce/InterfaceManagers/GhInterfaceManager.hpp"
#include "Evolution/Systems/Cce/InterfaceManagers/GhLocalTimeStepping.hpp"
#include "Evolution/Systems/Cce/InterfaceManagers/GhLockstep.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/System.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "Evolution/Systems/Cce/WorldtubeBufferUpdater.hpp"
#include "Evolution/Systems/Cce/WorldtubeDataManager.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Actions/NumericInitialData.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Bjorhus.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/DirichletMinkowski.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryConditions/Outflow.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/RegisterDerived.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/ConstraintDamping/RegisterDerivedWithCharm.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Equations.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/InitializeDampedHarmonic.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Initialize.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Evolution/TypeTraits.hpp"
#include "IO/Importers/Actions/RegisterWithElementDataReader.hpp"
#include "IO/Importers/ElementDataReader.hpp"
#include "IO/Observer/Actions/ObserverRegistration.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags.hpp"
#include "NumericalAlgorithms/Interpolation/BarycentricRationalSpanInterpolator.hpp"
#include "NumericalAlgorithms/Interpolation/CubicSpanInterpolator.hpp"
#include "NumericalAlgorithms/Interpolation/LinearSpanInterpolator.hpp"
#include "NumericalAlgorithms/Interpolation/SpanInterpolator.hpp"
#include "NumericalAlgorithms/LinearOperators/ExponentialFilter.hpp"
#include "NumericalAlgorithms/LinearOperators/FilterAction.hpp"
#include "Options/FactoryHelpers.hpp"
#include "Options/Options.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseControl/CheckpointAndExitAfterWallclock.hpp"
#include "Parallel/PhaseControl/ExecutePhaseChange.hpp"
#include "Parallel/PhaseControl/VisitAndReturn.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Actions/MemoryMonitor/ContributeMemoryData.hpp"
#include "ParallelAlgorithms/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "ParallelAlgorithms/Actions/SetupDataBox.hpp"
#include "ParallelAlgorithms/Actions/TerminatePhase.hpp"
#include "ParallelAlgorithms/Events/Factory.hpp"
#include "ParallelAlgorithms/Events/MonitorMemory.hpp"
#include "ParallelAlgorithms/Events/ObserveTimeStep.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/CleanUpInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/ElementInitInterpPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolationTargetReceiveVars.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceivePoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceiveVolumeData.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorRegisterElement.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/TryToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ErrorOnFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/IgnoreFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveSurfaceData.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/SendGhWorldtubeData.hpp"
#include "ParallelAlgorithms/Interpolation/Events/Interpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Events/InterpolateWithoutInterpComponent.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Interpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/ApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/Sphere.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/SphericalKerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/DetAndInverseSpatialMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ConstraintGammas.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "PointwiseFunctions/GeneralRelativity/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Actions/AdvanceTime.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Actions/RecordTimeStepperData.hpp"
#include "Time/Actions/SelfStartActions.hpp"
#include "Time/Actions/UpdateU.hpp"
#include "Time/StepChoosers/Cfl.hpp"
#include "Time/StepChoosers/Constant.hpp"
#include "Time/StepChoosers/ErrorControl.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/StepChoosers/Increase.hpp"
#include "Time/StepChoosers/PreventRapidIncrease.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/StepChoosers/StepToTimes.hpp"
#include "Time/StepControllers/Factory.hpp"
#include "Time/StepControllers/StepController.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeSequence.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
// IWYU pragma: no_forward_declare MathFunction
struct Inertial;
}  // namespace Frame
namespace PUP {
class er;
}  // namespace PUP
namespace Parallel {
template <typename Metavariables>
class CProxy_GlobalCache;
}  // namespace Parallel
/// \endcond

namespace detail {
template <typename InitialData>
constexpr auto make_default_phase_order() {
  if constexpr (evolution::is_numeric_initial_data_v<InitialData>) {
    return std::array{Parallel::Phase::Initialization,
                      Parallel::Phase::RegisterWithElementDataReader,
                      Parallel::Phase::ImportInitialData,
                      Parallel::Phase::InitializeInitialDataDependentQuantities,
                      Parallel::Phase::Register,
                      Parallel::Phase::InitializeTimeStepperHistory,
                      Parallel::Phase::Evolve,
                      Parallel::Phase::Exit};
  } else {
    return std::array{Parallel::Phase::Initialization,
                      Parallel::Phase::InitializeInitialDataDependentQuantities,
                      Parallel::Phase::Register,
                      Parallel::Phase::InitializeTimeStepperHistory,
                      Parallel::Phase::Evolve,
                      Parallel::Phase::Exit};
  }
}
}  // namespace detail

template <typename InitialData, typename BoundaryConditions>
struct EvolutionMetavars<3, InitialData, BoundaryConditions>
    : CharacteristicExtractDefaults {
  template <bool DuringSelfStart>
  struct CceWorldtubeTarget;

  static constexpr size_t volume_dim = 3;
  static constexpr bool use_damped_harmonic_rollon = true;
  using initial_data = InitialData;
  using system = GeneralizedHarmonic::System<volume_dim>;
  static constexpr dg::Formulation dg_formulation =
      dg::Formulation::StrongInertial;
  using temporal_id = Tags::TimeStepId;
  static constexpr bool local_time_stepping = false;
  // Set override_functions_of_time to true to override the
  // 2nd or 3rd order piecewise polynomial functions of time using
  // `read_spec_piecewise_polynomial()`
  static constexpr bool override_functions_of_time = false;

  using analytic_solution_fields = typename system::variables_tag::tags_list;

  using initialize_initial_data_dependent_quantities_actions =
      tmpl::list<GeneralizedHarmonic::gauges::Actions::InitializeDampedHarmonic<
                     volume_dim, use_damped_harmonic_rollon>,
                 Parallel::Actions::TerminatePhase>;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
  struct domain : tt::ConformsTo<::domain::protocols::Metavariables> {
    static constexpr bool enable_time_dependent_maps = true;
  };
  using analytic_solution_tag = Tags::AnalyticSolution<BoundaryConditions>;

  using analytic_compute =
      evolution::Tags::AnalyticSolutionsCompute<volume_dim,
                                                analytic_solution_fields>;
  using deriv_compute = ::Tags::DerivCompute<
      typename system::variables_tag,
      ::domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                      Frame::Inertial>,
      typename system::gradient_variables>;
  using error_compute = Tags::ErrorsCompute<analytic_solution_fields>;
  using error_tags = db::wrap_tags_in<Tags::Error, analytic_solution_fields>;

  struct AhA : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::Time;
    using tags_to_observe = ::ah::tags_for_observing<::Frame::Inertial>;
    using surface_tags_to_observe = ::ah::surface_tags_for_observing;
    using compute_vars_to_interpolate = ah::ComputeHorizonVolumeQuantities;
    using vars_to_interpolate_to_target =
        ::ah::vars_to_interpolate_to_target<volume_dim, ::Frame::Inertial>;
    using compute_items_on_target =
        ::ah::compute_items_on_target<volume_dim, ::Frame::Inertial>;
    using compute_target_points =
        intrp::TargetPoints::ApparentHorizon<AhA, ::Frame::Inertial>;
    using post_interpolation_callback =
        intrp::callbacks::FindApparentHorizon<AhA, ::Frame::Inertial>;
    using horizon_find_failure_callback =
        intrp::callbacks::IgnoreFailedApparentHorizon;
    using post_horizon_find_callbacks = tmpl::list<
        intrp::callbacks::ObserveTimeSeriesOnSurface<tags_to_observe, AhA>,
        intrp::callbacks::ObserveSurfaceData<surface_tags_to_observe, AhA>>;
  };

  using interpolator_source_vars = ::ah::source_vars<volume_dim>;

  using observe_fields = tmpl::append<
      tmpl::push_back<
          analytic_solution_fields, gr::Tags::Lapse<DataVector>,
          GeneralizedHarmonic::Tags::GaugeConstraintCompute<volume_dim,
                                                            ::Frame::Inertial>,
          GeneralizedHarmonic::Tags::TwoIndexConstraintCompute<
              volume_dim, ::Frame::Inertial>,
          GeneralizedHarmonic::Tags::ThreeIndexConstraintCompute<
              volume_dim, ::Frame::Inertial>,
          GeneralizedHarmonic::Tags::DerivSpatialMetricCompute<
              volume_dim, ::Frame::Inertial>,
          gr::Tags::SpatialChristoffelFirstKindCompute<
              volume_dim, ::Frame::Inertial, DataVector>,
          gr::Tags::SpatialChristoffelSecondKindCompute<
              volume_dim, ::Frame::Inertial, DataVector>,
          ::Tags::DerivTensorCompute<
              gr::Tags::SpatialChristoffelSecondKind<
                  volume_dim, ::Frame::Inertial, DataVector>,
              ::domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                              Frame::Inertial>>,
          gr::Tags::SpatialRicciCompute<volume_dim, ::Frame::Inertial,
                                        DataVector>,
          gr::Tags::SpatialRicciScalarCompute<volume_dim, ::Frame::Inertial,
                                              DataVector>,
          // following tags added to observe constraints
          ::Tags::PointwiseL2NormCompute<
              GeneralizedHarmonic::Tags::GaugeConstraint<volume_dim,
                                                         ::Frame::Inertial>>,
          ::Tags::PointwiseL2NormCompute<
              GeneralizedHarmonic::Tags::TwoIndexConstraint<volume_dim,
                                                            ::Frame::Inertial>>,
          ::Tags::PointwiseL2NormCompute<
              GeneralizedHarmonic::Tags::ThreeIndexConstraint<
                  volume_dim, ::Frame::Inertial>>,
          ::domain::Tags::Coordinates<volume_dim, Frame::Grid>,
          ::domain::Tags::Coordinates<volume_dim, Frame::Inertial>>,
      error_tags,
      // The 4-index constraint is only implemented in 3d
      tmpl::conditional_t<
          volume_dim == 3,
          tmpl::list<
              GeneralizedHarmonic::Tags::
                  FourIndexConstraintCompute<3, ::Frame::Inertial>,
              GeneralizedHarmonic::Tags::FConstraintCompute<3,
                                                            ::Frame::Inertial>,
              ::Tags::PointwiseL2NormCompute<
                  GeneralizedHarmonic::Tags::FConstraint<3, ::Frame::Inertial>>,
              ::Tags::PointwiseL2NormCompute<
                  GeneralizedHarmonic::Tags::FourIndexConstraint<
                      3, ::Frame::Inertial>>,
              GeneralizedHarmonic::Tags::ConstraintEnergyCompute<
                  3, ::Frame::Inertial>>,
          tmpl::list<>>>;
  using non_tensor_compute_tags =
      tmpl::list<::Events::Tags::ObserverMeshCompute<volume_dim>,
                 analytic_compute, error_compute,
                 ::Tags::TimeAndPreviousCompute>;

  using cce_vars_to_interpolate_to_target = interpolator_source_vars;

  // the subset of step choosers that are valid for the DG system only
  using dg_step_choosers = StepChoosers::standard_step_choosers<system>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<DenseTrigger, DenseTriggers::standard_dense_triggers>,
        tmpl::pair<DomainCreator<volume_dim>, domain_creators<volume_dim>>,
        tmpl::pair<
            Event,
            tmpl::flatten<tmpl::list<
                intrp::Events::Interpolate<3, AhA, interpolator_source_vars>,
                Events::Completion,
                intrp::Events::InterpolateWithoutInterpComponent<
                    volume_dim, CceWorldtubeTarget<true>, EvolutionMetavars,
                    cce_vars_to_interpolate_to_target>,
                intrp::Events::Interpolate<volume_dim,
                                           CceWorldtubeTarget<false>,
                                           cce_vars_to_interpolate_to_target>,
                dg::Events::field_observations<volume_dim, Tags::Time,
                                               observe_fields,
                                               non_tensor_compute_tags>,
                Events::time_events<system>>>>,
        tmpl::pair<GeneralizedHarmonic::BoundaryConditions::BoundaryCondition<
                       volume_dim>,
                   GeneralizedHarmonic::BoundaryConditions::
                       standard_boundary_conditions<volume_dim>>,
        tmpl::pair<LtsTimeStepper, TimeSteppers::lts_time_steppers>,
        tmpl::pair<
            PhaseChange,
            tmpl::list<
                PhaseControl::VisitAndReturn<Parallel::Phase::LoadBalancing>,
                PhaseControl::VisitAndReturn<Parallel::Phase::WriteCheckpoint>,
                PhaseControl::CheckpointAndExitAfterWallclock>>,
        tmpl::pair<StepChooser<StepChooserUse::LtsStep>,
                   tmpl::remove_duplicates<tmpl::append<
                       StepChoosers::standard_step_choosers<system>,
                       cce_step_choosers>>>,
        tmpl::pair<
            StepChooser<StepChooserUse::Slab>,
            StepChoosers::standard_slab_choosers<system, local_time_stepping>>,
        tmpl::pair<StepController, StepControllers::standard_step_controllers>,
        tmpl::pair<TimeSequence<double>,
                   TimeSequences::all_time_sequences<double>>,
        tmpl::pair<TimeSequence<std::uint64_t>,
                   TimeSequences::all_time_sequences<std::uint64_t>>,
        tmpl::pair<TimeStepper, TimeSteppers::time_steppers>,
        tmpl::pair<Trigger, tmpl::append<Triggers::logical_triggers,
                                         Triggers::time_triggers>>>;
  };

  // A tmpl::list of tags to be added to the GlobalCache by the
  // metavariables
  using const_global_cache_tags = tmpl::list<
      analytic_solution_tag,
      GeneralizedHarmonic::ConstraintDamping::Tags::DampingFunctionGamma0<
          volume_dim, Frame::Grid>,
      GeneralizedHarmonic::ConstraintDamping::Tags::DampingFunctionGamma1<
          volume_dim, Frame::Grid>,
      GeneralizedHarmonic::ConstraintDamping::Tags::DampingFunctionGamma2<
          volume_dim, Frame::Grid>>;

  using dg_registration_list =
      tmpl::list<observers::Actions::RegisterEventsWithObservers,
                 intrp::Actions::RegisterElementWithInterpolator>;

  static constexpr auto default_phase_order =
      detail::make_default_phase_order<initial_data>();

  template <bool DuringSelfStart>
  using step_actions = tmpl::list<
      tmpl::conditional_t<DuringSelfStart,
                          tmpl::list<Cce::Actions::InterpolateDuringSelfStart<
                              CceWorldtubeTarget<true>>>,
                          tmpl::list<>>,
      evolution::dg::Actions::ComputeTimeDerivative<EvolutionMetavars>,
      tmpl::conditional_t<
          local_time_stepping,
          tmpl::list<evolution::Actions::RunEventsAndDenseTriggers<>,
                     evolution::dg::Actions::ApplyLtsBoundaryCorrections<
                         EvolutionMetavars>>,
          tmpl::list<
              evolution::dg::Actions::ApplyBoundaryCorrectionsToTimeDerivative<
                  EvolutionMetavars>,
              Actions::RecordTimeStepperData<>,
              evolution::Actions::RunEventsAndDenseTriggers<>,
              Actions::UpdateU<>,
              dg::Actions::Filter<
                  Filters::Exponential<0>,
                  tmpl::list<gr::Tags::SpacetimeMetric<
                                 volume_dim, Frame::Inertial, DataVector>,
                             GeneralizedHarmonic::Tags::Pi<volume_dim,
                                                           Frame::Inertial>,
                             GeneralizedHarmonic::Tags::Phi<
                                 volume_dim, Frame::Inertial>>>>>>;

  using initialization_actions = tmpl::list<
      Actions::SetupDataBox,
      Initialization::Actions::TimeAndTimeStep<EvolutionMetavars>,
      evolution::dg::Initialization::Domain<volume_dim,
                                            override_functions_of_time>,
      Initialization::Actions::NonconservativeSystem<system>,
      std::conditional_t<
          evolution::is_numeric_initial_data_v<initial_data>, tmpl::list<>,
          evolution::Initialization::Actions::SetVariables<
              ::domain::Tags::Coordinates<volume_dim, Frame::ElementLogical>>>,
      Initialization::Actions::AddComputeTags<::Tags::DerivCompute<
          typename system::variables_tag,
          ::domain::Tags::InverseJacobian<volume_dim, Frame::ElementLogical,
                                          Frame::Inertial>,
          typename system::gradient_variables>>,
      Initialization::Actions::TimeStepperHistory<EvolutionMetavars>,
      GeneralizedHarmonic::Actions::InitializeGhAnd3Plus1Variables<volume_dim>,
      Initialization::Actions::AddComputeTags<tmpl::push_back<
          StepChoosers::step_chooser_compute_tags<EvolutionMetavars>>>,
      ::evolution::dg::Initialization::Mortars<volume_dim, system>,
      evolution::Actions::InitializeRunEventsAndDenseTriggers,
      intrp::Actions::ElementInitInterpPoints<
          intrp::Tags::InterpPointInfo<EvolutionMetavars>>,
      Initialization::Actions::RemoveOptionsAndTerminatePhase>;

  using gh_dg_element_array = DgElementArray<
      EvolutionMetavars,
      tmpl::flatten<tmpl::list<
          Parallel::PhaseActions<Parallel::Phase::Initialization,
                                 initialization_actions>,
          tmpl::conditional_t<
              evolution::is_numeric_initial_data_v<initial_data>,
              tmpl::list<
                  Parallel::PhaseActions<
                      Parallel::Phase::RegisterWithElementDataReader,
                      tmpl::list<
                          importers::Actions::RegisterWithElementDataReader,
                          Parallel::Actions::TerminatePhase>>,
                  Parallel::PhaseActions<
                      Parallel::Phase::ImportInitialData,
                      tmpl::list<
                          GeneralizedHarmonic::Actions::ReadNumericInitialData<
                              evolution::OptionTags::NumericInitialData>,
                          GeneralizedHarmonic::Actions::SetNumericInitialData<
                              evolution::OptionTags::NumericInitialData>,
                          Parallel::Actions::TerminatePhase>>>,
              tmpl::list<>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeInitialDataDependentQuantities,
              initialize_initial_data_dependent_quantities_actions>,
          Parallel::PhaseActions<Parallel::Phase::Register,
                                 tmpl::list<dg_registration_list,
                                            Parallel::Actions::TerminatePhase>>,
          Parallel::PhaseActions<
              Parallel::Phase::InitializeTimeStepperHistory,
              SelfStart::self_start_procedure<step_actions<true>, system>>,
          Parallel::PhaseActions<
              Parallel::Phase::Evolve,
              tmpl::list<Actions::RunEventsAndTriggers, Actions::ChangeSlabSize,
                         step_actions<false>, Actions::AdvanceTime,
                         PhaseControl::Actions::ExecutePhaseChange>>>>>;

  template <bool DuringSelfStart>
  struct CceWorldtubeTarget
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = tmpl::conditional_t<DuringSelfStart, ::Tags::TimeStepId,
                                            ::Tags::TimeAndPrevious>;

    static std::string name() {
      return DuringSelfStart ? "SelfStartCceWorldtubeTarget"
                             : "CceWorldtubeTarget";
    }
    using compute_items_on_source = tmpl::list<>;
    using compute_items_on_target = tmpl::list<>;
    using compute_target_points =
        intrp::TargetPoints::Sphere<CceWorldtubeTarget, ::Frame::Inertial>;
    using post_interpolation_callback = intrp::callbacks::SendGhWorldtubeData<
        Cce::CharacteristicEvolution<EvolutionMetavars>, DuringSelfStart>;
    using vars_to_interpolate_to_target = cce_vars_to_interpolate_to_target;
    template <typename Metavariables>
    using interpolating_component = gh_dg_element_array;
  };

  using interpolation_target_tags =
      tmpl::list<AhA, CceWorldtubeTarget<true>, CceWorldtubeTarget<false>>;

  using observed_reduction_data_tags = observers::collect_reduction_data_tags<
      tmpl::at<typename factory_creation::factory_classes, Event>>;

  using cce_boundary_component = Cce::GhWorldtubeBoundary<EvolutionMetavars>;

  template <typename ParallelComponent>
  struct registration_list {
    using type = std::conditional_t<
        std::is_same_v<ParallelComponent, gh_dg_element_array>,
        dg_registration_list, tmpl::list<>>;
  };

  using component_list = tmpl::flatten<tmpl::list<
      observers::Observer<EvolutionMetavars>,
      observers::ObserverWriter<EvolutionMetavars>,
      std::conditional_t<evolution::is_numeric_initial_data_v<InitialData>,
                         importers::ElementDataReader<EvolutionMetavars>,
                         tmpl::list<>>,
      intrp::Interpolator<EvolutionMetavars>,
      tmpl::transform<interpolation_target_tags,
                      tmpl::bind<intrp::InterpolationTarget,
                                 tmpl::pin<EvolutionMetavars>, tmpl::_1>>,
      cce_boundary_component, gh_dg_element_array,
      Cce::CharacteristicEvolution<EvolutionMetavars>>>;

  static constexpr Options::String help{"GH Single BH + CCE"};
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling,
    &setup_memory_allocation_failure_reporting,
    &disable_openblas_multithreading,
    &::domain::creators::time_dependence::register_derived_with_charm,
    &::domain::FunctionsOfTime::register_derived_with_charm,
    &GeneralizedHarmonic::BoundaryCorrections::register_derived_with_charm,
    &Cce::register_initialize_j_with_charm<
        metavariables::uses_partially_flat_cartesian_coordinates>,
    &Parallel::register_derived_classes_with_charm<Cce::WorldtubeDataManager>,
    &Parallel::register_derived_classes_with_charm<intrp::SpanInterpolator>,
    &::domain::creators::register_derived_with_charm,
    &GeneralizedHarmonic::ConstraintDamping::register_derived_with_charm,
    &Parallel::register_factory_classes_with_charm<metavariables>};

static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
