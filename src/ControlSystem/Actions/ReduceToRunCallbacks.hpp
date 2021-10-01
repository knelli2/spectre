// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <ostream>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/RunCallbacks.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Printf.hpp"
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace control_system {
struct SubTrackTranslation;
}  // namespace control_system
/// \endcond

namespace control_system {
namespace Actions {
template <typename ControlSystems>
struct ReduceToRunCallbacks {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(const db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const LinkedMessageId<double>& measurement_id,
                    const std::vector<double>& integrals,
                    const std::vector<DataVector>& vec_coords) noexcept {
    // const typename ::Tags::MeasureTranslationResult::type&
    //    data_from_element) noexcept {
    std::ostringstream os;
    // os << "Reduction finished! Inside ReduceToRunCallbacks.\n";
    // Parallel::printf(os.str());
    double coord_min = std::numeric_limits<double>::max();
    double coord_max = std::numeric_limits<double>::min();
    for (auto& dv : vec_coords) {
      if (min(dv) < coord_min) {
        coord_min = min(dv);
      }
      if (max(dv) > coord_max) {
        coord_max = max(dv);
      }
    }
    const double center_coord = 0.5 * (coord_max + coord_min);
    std::vector<double> integrals2 = integrals;
    integrals2.push_back(center_coord);
    const auto temp_box = db::create<
        db::AddSimpleTags<control_system::Tags::MeasureTranslationResult>>(
        integrals2);
    //        data_from_element);
    control_system::RunCallbacks<control_system::SubTrackTranslation,
                                 ControlSystems>::apply(temp_box, cache,
                                                        measurement_id);
  }
};
}  // namespace Actions
}  // namespace control_system
