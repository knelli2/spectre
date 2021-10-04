// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ControlSystem/ControlErrors/DataVectorHelpers.hpp"
#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim>
struct Domain;
}  // namespace domain::Tags
/// \endcond

namespace {
struct SomeTagHere;
}

namespace control_system {
namespace ControlErrors {
struct Rotation {
  template <size_t DerivOrder, typename Metavariables,
            typename... TupleTags>
  static DataVector apply(
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& function_of_time_name,
      const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& domain = get<domain::Tags::Domain<DerivOrder - 1>>(cache);

    const double grid_position_of_A = domain.x_coord_object_A();
    const double grid_position_of_B = domain.x_coord_object_B();
    const DataVector grid_separation{
        {grid_position_of_A - grid_position_of_B, 0.0, 0.0}};
    const DataVector& current_position_of_A = get<SomeTagHere>(measurements);
    const DataVector& current_position_of_B = get<SomeTagHere>(measurements);
    const DataVector current_separation =
        current_position_of_A - current_position_of_B;

    const double current_dot_grid_separation =
        dot(current_separation, grid_separation);
    const DataVector current_cross_grid_separation =
        cross(current_separation, grid_separation);

    return -current_cross_grid_sepatation / current_dot_grid_separation;
  }
};
}  // namespace ControlErrors
}  // namespace control_system
