// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

#include "ControlSystem/Component.hpp"
#include "ControlSystem/Protocols/Measurement.hpp"
#include "ControlSystem/Protocols/Submeasurement.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/Actions/FunctionsOfTimeAreReady.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace control_system::Tags {
struct WriteDataToDisk;
}  // namespace control_system::Tags
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace Tags {
struct PreviousTriggerTime;
}  // namespace Tags
namespace ScalarWave::Tags {
struct Psi;
}  // namespace ScalarWave::Tags
/// \endcond

namespace control_system {
namespace measurements {
/// \cond
template <typename ControlSystems>
struct PostReductionSendWaveRadiusToControlSystem;
/// \endcond
}  // namespace measurements

/*!
 * \brief An `::Event` that computes the center of mass for $x > 0$ and $x < 0$
 * where $x$ is the in the `Frame::Grid`.
 *
 * \details See
 * `control_system::measurements::center_of_mass_integral_on_element` for the
 * calculation of the CoM. This event then does a reduction and calls
 * `control_system::PostReductionSendWaveRadiusToControlSystem` as a post
 * reduction callback.
 *
 * \tparam ControlSystems `tmpl::list` of all control systems that use this
 * event.
 */
template <typename ControlSystems>
class ScalarWaveTracker : public ::Event {
 public:
  /// \cond
  // LCOV_EXCL_START
  explicit ScalarWaveTracker(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ScalarWaveTracker);  // NOLINT
  // LCOV_EXCL_STOP
  /// \endcond

  // This event is created during control system initialization, not
  // from the input file.
  static constexpr bool factory_creatable = false;
  ScalarWaveTracker() = default;

  using compute_tags_for_observation_box = tmpl::list<>;

  using return_tags = tmpl::list<>;
  using argument_tags =
      tmpl::list<::Tags::Time, ::Tags::PreviousTriggerTime,
                 ScalarWave::Tags::Psi,
                 ::Events::Tags::ObserverCoordinates<3, Frame::Grid>>;

  template <typename Metavariables, typename ArrayIndex,
            typename ParallelComponent>
  void operator()(const double time, const std::optional<double>& previous_time,
                  const Scalar<DataVector>& psi,
                  const tnsr::I<DataVector, 3, Frame::Grid>& x_grid,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ArrayIndex& array_index,
                  const ParallelComponent* const /*meta*/,
                  const ObservationValue& /*observation_value*/) const {
    const LinkedMessageId<double> measurement_id{time, previous_time};
    const Scalar<DataVector> grid_radius = magnitude(x_grid);
    const auto iterator_to_max_element = alg::max_element(get(psi));
    const auto index_of_max_element = static_cast<size_t>(
        std::distance(get(psi).begin(), iterator_to_max_element));

    ASSERT(index_of_max_element < get(grid_radius).size(),
           "Index of max element " << index_of_max_element
                                   << " outside the range of the coords vector "
                                   << get(grid_radius).size() - 1);

    // Reduction
    auto my_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache)[array_index];
    // We need a place to run RunCallback on... this does not need to be
    // the control system using the CoM data.
    auto& reduction_target_proxy = Parallel::get_parallel_component<
        ControlComponent<Metavariables, tmpl::front<ControlSystems>>>(cache);
    Parallel::ReductionData<
        Parallel::ReductionDatum<LinkedMessageId<double>, funcl::AssertEqual<>>,
        Parallel::ReductionDatum<double, funcl::Max<>>>
        reduction_data{measurement_id, get(grid_radius)[index_of_max_element]};
    Parallel::contribute_to_reduction<
        measurements::PostReductionSendWaveRadiusToControlSystem<
            ControlSystems>>(std::move(reduction_data), my_proxy,
                             reduction_target_proxy);
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*component*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }
};

/// \cond
template <typename ControlSystems>
PUP::able::PUP_ID ScalarWaveTracker<ControlSystems>::my_PUP_ID = 0;  // NOLINT
/// \endcond

namespace measurements {
namespace Tags {
struct WaveRadius : db::SimpleTag {
  using type = double;
};
}  // namespace Tags

struct ScalarWavePeakRadius : tt::ConformsTo<protocols::Measurement> {
  struct PeakRadius : tt::ConformsTo<protocols::Submeasurement> {
    static std::string name() { return "ScalarWavePeakRadius::PeakRadius"; }
    /// Unused tag needed to conform to the submeasurement protocol.
    template <typename ControlSystems>
    using interpolation_target_tag = void;

    template <typename ControlSystems>
    using event = ScalarWaveTracker<ControlSystems>;
  };
  /// List of submeasurements used by this measurement -- only WaveRadius
  /// here.
  using submeasurements = tmpl::list<PeakRadius>;
};

template <typename ControlSystems>
struct PostReductionSendWaveRadiusToControlSystem {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const LinkedMessageId<double>& id, const double radius) {
    const auto center_databox =
        db::create<db::AddSimpleTags<Tags::WaveRadius>>(radius);
    // Send results to the control system(s)
    RunCallbacks<ScalarWavePeakRadius::PeakRadius, ControlSystems>::apply(
        center_databox, cache, id);

    if (Parallel::get<control_system::Tags::WriteDataToDisk>(cache)) {
      std::vector<double> data_to_write{id.id, radius};

      auto& writer_proxy = Parallel::get_parallel_component<
          observers::ObserverWriter<Metavariables>>(cache);

      Parallel::threaded_action<
          observers::ThreadedActions::WriteReductionDataRow>(
          // Node 0 is always the writer
          writer_proxy[0], subfile_path_, legend_,
          std::make_tuple(std::move(data_to_write)));
    }
  }

 private:
  const static inline std::vector<std::string> legend_{"Time",
                                                       "WavePeakRadius"};
  const static inline std::string subfile_path_{
      "/ControlSystems/ScalarWaveRadius"};
};
}  // namespace measurements
}  // namespace control_system
