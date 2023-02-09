// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "ApparentHorizons/ComputeHorizonVolumeQuantities.hpp"
#include "ApparentHorizons/HorizonAliases.hpp"
#include "ApparentHorizons/ObserveCenters.hpp"
#include "ApparentHorizons/Tags.hpp"
#include "ControlSystem/Protocols/Measurement.hpp"
#include "ControlSystem/Protocols/Submeasurement.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/LinkedMessageId.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/ObjectLabel.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ErrorOnFailedApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/FindApparentHorizon.hpp"
#include "ParallelAlgorithms/Interpolation/Interpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/ApparentHorizon.hpp"
#include "Time/Tags.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <size_t VolumeDim>
class ElementId;
template <size_t Dim>
class Mesh;
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace domain::Tags {
template <size_t Dim>
struct Mesh;
}  // namespace domain::Tags
/// \endcond

namespace control_system {
struct BothNSCenters : tt::ConformsTo<protocols::Measurement> {
  template <::domain::ObjectLabel Horizon>
  struct FindCenter : tt::ConformsTo<protocols::Submeasurement> {
    using argument_tags = tmpl::list<list_of_tags_needed_to_do_integrations>;

    template <typename Metavariables, typename ParallelComponent,
              typename ControlSystems>
    static void apply(
        const Mesh<3>& mesh,
        const tnsr::aa<DataVector, 3, ::Frame::Inertial>& tag_types...,
        const LinkedMessageId<double>& measurement_id,
        Parallel::GlobalCache<Metavariables>& cache,
        const ElementId<3>& array_index,
        const ParallelComponent* const /*meta*/, ControlSystems /*meta*/) {
      // This is run on the elements.
      // 1. Do integration
      // 2. call
      //     Parallel::reduction<ActionToBeCalledAfterReduction>(
      //      proxy_of_control_system, cache, index,
      //      integrated_values_neg_x, integrated_values_pos_x);
    }
  };

  using submeasurements = tmpl::list<FindCenter<::domain::ObjectLabel::A>,
                                     FindCenter<::domain::ObjectLabel::B>>;
};

struct ActionToBeCalledAfterReduction {
  static void apply(cache, index, total_integrated_value_neg_x,
                    integrated_value_pos_x) {
    // This is run on proxy defined above
    // total_integrated_value is the entire integral
    // 1. Do math
    // create temporary box with just necessary tags???
    // db::create<>
    // 2. RunCallbacks::apply(box, cache, centers)
  }
};
}  // namespace control_system
