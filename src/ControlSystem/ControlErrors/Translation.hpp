// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/math/quaternion.hpp>

#include "ControlSystem/ControlErrors/DataVectorHelpers.hpp"
#include "ControlSystem/KyleExample/SubTrackTranslation.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/QuaternionHelpers.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system {
namespace ControlErrors {
struct Translation {
  template <size_t DerivOrder, typename Metavariables, typename... TupleTags>
  static DataVector apply(
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& function_of_time_name,
      const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& domain = get<domain::Tags::Domain<DerivOrder - 1>>(cache);
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);

    using quat = boost::math::quaternion<double>;

    const quat quaternion = datavector_to_quaternion(
        functions_of_time.at("Rotation")->func(time)[0]);
    const double expansion_coef =
        functions_of_time.at("Expansion")->func(time)[0][0];

    const DataVector grid_position_of_A = {
        {domain.x_coord_object_A(), 0.0, 0.0}};
    const DataVector grid_position_of_B = {
        {domain.x_coord_object_B(), 0.0, 0.0}};
    const DataVector grid_separation = grid_position_of_A - grid_position_of_B;
    const DataVector& current_position_of_A = get<SomeTagHere>(measurements);
    const DataVector& current_position_of_B = get<SomeTagHere>(measurements);
    const DataVector current_separation =
        current_position_of_A - current_position_of_B;

    const double current_dot_grid_separation =
        dot(current_separation, grid_separation);
    const DataVector current_cross_grid_separation =
        cross(current_separation, grid_separation);

    const DataVector rotation_param = datavector_to_quaternion(
        cross(current_cross_grid_separation / current_dot_grid_separation,
              grid_position_of_A));

    const double expansion_param =
        expansion_coef *
        ((current_position_of_A[0] - current_position_of_B[0]) /
             (grid_position_of_A[0] - grid_position_of_B[0]) -
         1.0);

    // From eq 42 in 1304.3067
    const quat temp_middle = datavector_to_quaternion(
        current_position_of_A -
        (1.0 + expansion_param / expansion_coef) * grid_position_of_A -
        rotation_param);
    return expansion_coef * quaternion_to_datavector(quaternion * temp_middle *
                                                     conj(quaternion));
  }
};

struct Translation1D {
  template <size_t DerivOrder, typename... TupleTags>
  static DataVector apply(
      const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& integral_vector =
        get<control_system::QueueTags::MeasureTranslation>(measurements);
    DataVector target_value{DerivOrder - 1, 0.0};
    // integral(psi^2 * x_i) / integral(psi^2)
    // basically a CoM
    for (size_t i = 1; i < DerivOrder; i++) {
      target_value[i - 1] = integral_vector[i] / integral_vector[0];
    }

    const double& center_of_interval = integral_vector[DerivOrder];
    return target_value - center_of_interval;
  }
};
}  // namespace ControlErrors
}  // namespace control_system
