// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/math/quaternion.hpp>
#include <cstddef>
#include <pup.h>

#include "ControlSystem/ControlErrors/Expansion.hpp"
#include "ControlSystem/ControlErrors/Rotation.hpp"
#include "ControlSystem/DataVectorHelpers.hpp"
#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/QuaternionHelpers.hpp"
#include "Options/Options.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace control_system::Tags {
struct MeasurementTimescales;
}  // namespace control_system::Tags
/// \endcond

namespace control_system {
namespace ControlErrors {
/*!
 * \brief Control error in the 3D \link
 * domain::CoordinateMaps::TimeDependent::Translation Translation \endlink
 * coordinate map
 *
 * \details Computes the error in how much the system has translated by using
 * Eq. (42) from \cite Ossokine2013zga. The equation is
 *
 * \f[ \left(0, \delta\vec{T}\right) = a\mathbf{q}\left(\mathbf{x}_B -
 * \mathbf{c}_B - \mathbf{\delta q}\wedge\mathbf{c}_B - \frac{\delta
 * a}{a}\mathbf{c}_B \right)\mathbf{q}^*
 * \f]
 *
 * where object B is located on the positive x-axis, bold face letters are
 * quaternions, vectors are promoted to quaternions as \f$ \mathbf{v} = (0,
 * \vec{v}) \f$, \f$ \mathbf{q} \f$ is the quaternion from the \link
 * domain::CoordinateMaps::TimeDependent::Rotation Rotation \endlink map, \f$ a
 * \f$ is the function \f$ a(t) \f$ from the \link
 * domain::CoordinateMaps::TimeDependent::CubicScale CubicScale \endlink map,
 * \f$ \mathbf{\delta q}\wedge\mathbf{c}_B \equiv (0, \delta\vec{q} \times
 * \vec{c}_B) \f$, \f$ \delta\vec{q} \f$ is the \link
 * control_system::ControlErrors::Rotation Rotation \endlink control error, and
 * \f$ \delta a\f$ is the \link control_system::ControlErrors::Expansion
 * Expansion \endlink control error. It is assumed that the positions of the
 * excision spheres are exactly along the x-axis.
 *
 * Requirements:
 * - This control error requires that there be exactly two objects in the
 *   simulation
 * - Currently both these objects must be black holes
 * - Currently this control system can only be used with the \link
 *   control_system::Systems::Translation Translation \endlink control system
 */
struct Translation : tt::ConformsTo<protocols::ControlError> {
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Computes the control error for translation control. This should not "
      "take any options."};

  void pup(PUP::er& /*p*/) {}

  template <typename Metavariables, typename... TupleTags>
  DataVector operator()(const Parallel::GlobalCache<Metavariables>& cache,
                        const double time,
                        const std::string& function_of_time_name,
                        const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& domain = get<domain::Tags::Domain<3>>(cache);
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);

    using quat = boost::math::quaternion<double>;

    const quat quaternion = datavector_to_quaternion(
        functions_of_time.at("Rotation")->func(time)[0]);
    const double expansion_factor =
        functions_of_time.at("Expansion")->func(time)[0][0];

    using horizon_label = control_system::ah::HorizonLabel;
    using center_A = control_system::QueueTags::Center<horizon_label::AhA>;
    using center_B = control_system::QueueTags::Center<horizon_label::AhB>;

    const DataVector grid_position_of_B = array_to_datavector(
        domain.excision_spheres().at("ObjectBExcisionSphere").center());
    const DataVector current_position_of_B = get<center_B>(measurements);
    const DataVector grid_position_of_A = array_to_datavector(
        domain.excision_spheres().at("ObjectAExcisionSphere").center());
    const DataVector current_position_of_A = get<center_A>(measurements);

    const DataVector rotation_error =
        rotation_control_error_(cache, time, "Rotation", measurements);
    // Use B because it's on the positive x-axis, however, A would work as well.
    // Just so long as we are consistent.
    const DataVector rotation_error_cross_grid_pos_B =
        cross(rotation_error, grid_position_of_B);

    const double expansion_error =
        expansion_control_error_(cache, time, "Expansion", measurements)[0];

    // From eq. 42 in 1304.3067
    const quat middle_expression = datavector_to_quaternion(
        current_position_of_B -
        (1.0 + expansion_error / expansion_factor) * grid_position_of_B -
        rotation_error_cross_grid_pos_B);

    // Because we are converting from a quaternion to a DataVector, there will
    // be four components in the DataVector. However, rotation control only
    // requires three (the latter three to be exact, because the first component
    // should be 0. We ASSERT this also.)
    const DataVector result_with_four_components =
        expansion_factor *
        quaternion_to_datavector(quaternion * middle_expression *
                                 conj(quaternion));
    ASSERT(equal_within_roundoff(result_with_four_components[0], 0.0),
           "Error in computing translation control error. First component of "
           "resulting quaternion should be 0.0, but is " +
               get_output(result_with_four_components[0]) + " instead.");

    const auto trans_func_and_2_derivs =
        functions_of_time.at(function_of_time_name)->func_and_2_derivs(time);

    const auto& timescales =
        get<control_system::Tags::MeasurementTimescales>(cache);

    const auto timescale = timescales.at(function_of_time_name)->func(time)[0];

    Parallel::printf(
        "Time: %.17f\n"
        " grid_A: %s\n"
        " grid_B: %s\n"
        " current_A: %s\n"
        " current_B: %s\n"
        " rot_error: %s\n"
        " exp_error: %.17f\n"
        " trans_error: %s\n"
        " trans_func: %s\n"
        " trans_deriv: %s\n"
        " trans_2deriv: %s\n"
        " measure timescale: %s\n\n",
        time, grid_position_of_A, grid_position_of_B, current_position_of_A,
        current_position_of_B, rotation_error, expansion_error,
        result_with_four_components, trans_func_and_2_derivs[0],
        trans_func_and_2_derivs[1], trans_func_and_2_derivs[2], timescale);

    return {result_with_four_components[1], result_with_four_components[2],
            result_with_four_components[3]};
  }

 private:
  Rotation rotation_control_error_{};
  Expansion expansion_control_error_{};
};
}  // namespace ControlErrors
}  // namespace control_system
