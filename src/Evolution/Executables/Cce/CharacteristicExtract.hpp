// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Executables/Cce/CharacteristicExtractBase.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/BouncingBlackHole.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/GaugeWave.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/LinearizedBondiSachs.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/RobinsonTrautman.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/RotatingSchwarzschild.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/SphericalMetricData.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/TeukolskyWave.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/WorldtubeData.hpp"
#include "Evolution/Systems/Cce/BoundaryData.hpp"
#include "Evolution/Systems/Cce/Components/CharacteristicEvolution.hpp"
#include "Evolution/Systems/Cce/Components/WorldtubeBoundary.hpp"
#include "Evolution/Systems/Cce/Events/ObserveFields.hpp"
#include "Evolution/Systems/Cce/Events/ObserveTimeStep.hpp"
#include "Evolution/Systems/Cce/Initialize/ConformalFactor.hpp"
#include "Evolution/Systems/Cce/Initialize/InitializeJ.hpp"
#include "Evolution/Systems/Cce/Initialize/InverseCubic.hpp"
#include "Evolution/Systems/Cce/Initialize/NoIncomingRadiation.hpp"
#include "Evolution/Systems/Cce/Initialize/ZeroNonSmooth.hpp"
#include "Evolution/Systems/Cce/IntegrandInputSteps.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/System.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "Evolution/Systems/Cce/WorldtubeDataManager.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Options/String.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/Phase.hpp"
#include "ParallelAlgorithms/Events/Factory.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/EventsAndTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/LogicalTriggers.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/StepChoosers/Factory.hpp"
#include "Time/TimeSteppers/Factory.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Time/Triggers/TimeTriggers.hpp"
#include "Utilities/ProtocolHelpers.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

template <template <typename> class BoundaryComponent>
struct EvolutionMetavars : CharacteristicExtractDefaults<false> {
  using system = Cce::System<evolve_ccm>;
  static constexpr bool local_time_stepping = true;
  using cce_boundary_component = BoundaryComponent<EvolutionMetavars>;

  using component_list =
      tmpl::list<observers::ObserverWriter<EvolutionMetavars>,
                 cce_boundary_component,
                 Cce::CharacteristicEvolution<EvolutionMetavars>>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<LtsTimeStepper, TimeSteppers::lts_time_steppers>,
        tmpl::pair<StepChooser<StepChooserUse::LtsStep>, cce_step_choosers>,
        tmpl::pair<StepChooser<StepChooserUse::Slab>,
                   StepChoosers::standard_slab_choosers<
                       system, local_time_stepping, false>>,
        tmpl::pair<TimeSequence<double>,
                   TimeSequences::all_time_sequences<double>>,
        tmpl::pair<TimeSequence<std::uint64_t>,
                   TimeSequences::all_time_sequences<std::uint64_t>>,
        tmpl::pair<Event, tmpl::list<Cce::Events::ObserveFields,
                                     Cce::Events::ObserveTimeStep>>,
        tmpl::pair<Trigger, tmpl::append<Triggers::logical_triggers,
                                         Triggers::time_triggers>>>;
  };

  using observed_reduction_data_tags = tmpl::list<>;

  inline const static std::string help{
      "Perform Cauchy Characteristic Extraction using .h5 input data.\n"
      "Uses regularity-preserving formulation."};

  static constexpr std::array<Parallel::Phase, 4> default_phase_order{
      {Parallel::Phase::Initialization,
       Parallel::Phase::InitializeTimeStepperHistory, Parallel::Phase::Evolve,
       Parallel::Phase::Exit}};

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
