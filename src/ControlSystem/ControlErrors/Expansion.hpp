// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
/// \endcond

namespace {
struct SomeTagHere;
}

namespace control_system {
namespace ControlErrors {
struct Expansion {
  template <size_t DerivOrder, typename Metavariables,
            typename... TupleTags>
  static DataVector apply(
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& function_of_time_name,
      const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& domain = get<domain::Tags::Domain<DerivOrder-1>>(cache);
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);

    const double current_expansion_factor =
        functions_of_time.at(function_of_time_name)->func(time)[0][0];

    const double grid_position_of_A = domain.x_coord_object_A();
    const double grid_position_of_B = domain.x_coord_object_B();
    const double current_position_of_A = get<SomeTagHere>(measurements);
    const double current_position_of_B = get<SomeTagHere>(measurements);

    const double expected_expansion_factor =
        current_expansion_factor *
        (current_position_of_A - current_position_of_B) /
        (grid_position_of_A - grid_position_of_B);

    return expected_expansion_factor - current_expansion_factor;
  }
};
}  // namespace ControlErrors
}  // namespace control_system
