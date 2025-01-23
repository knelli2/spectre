// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>

#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Tags/QueueTags.hpp"
#include "ControlSystem/Tags/SystemTags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "Domain/Block.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/Tags/ObjectCenter.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
/// \endcond

namespace control_system::ControlErrors {
/*!
 * \brief Needs some dox
 */
struct InertialCenters : tt::ConformsTo<protocols::ControlError> {
  using object_centers =
      domain::object_list<domain::ObjectLabel::A, domain::ObjectLabel::B>;

  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Computes the control error for the inertial centers of two objects. "
      "This should not take any options."};

  void pup(PUP::er& /*p*/) {}

  template <typename Metavariables, typename... TupleTags>
  DataVector operator()(const ::TimescaleTuner<true>& /*unused*/,
                        const Parallel::GlobalCache<Metavariables>& cache,
                        const double time,
                        const std::string& function_of_time_name,
                        const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
    const auto& function_of_time = functions_of_time.at(function_of_time_name);
    const DataVector fot_inertial_positions = function_of_time->func(time)[0];

    const DataVector current_centers =
        functions_of_time.at(function_of_time_name)->func(time)[0];

    using center_A = control_system::QueueTags::Center<::domain::ObjectLabel::A,
                                                       Frame::Inertial>;
    using center_B = control_system::QueueTags::Center<::domain::ObjectLabel::B,
                                                       Frame::Inertial>;

    const auto& current_position_of_A = get<center_A>(measurements);
    const auto& current_position_of_B = get<center_B>(measurements);
    tnsr::I<DataVector, 3, Frame::Grid> current_positions{2_st, 0.0};
    for (size_t i = 0; i < 3; i++) {
      current_positions.get(i)[0] = current_position_of_A[i];
      current_positions.get(i)[1] = current_position_of_B[i];
    }

    const Domain<3>& domain = Parallel::get<domain::Tags::Domain<3>>(cache);
    const auto block_logical_coords =
        block_logical_coordinates(domain, current_positions);
    ASSERT(alg::all_of(block_logical_coords,
                       [](const auto& coord) { return coord.has_value(); }),
           "Centers are no longer in the domain.");

    DataVector control_error{6, 0.0};

    // Need to compare inertial position from the FoT to the inertial position
    // gotten by mapping the measured grid centers to the inertial frame
    for (size_t i = 0; i < 2; i++) {
      const tnsr::I<double, 3, Frame::Grid> current_position{std::array{
          i == 0 ? current_position_of_A[0] : current_position_of_B[0],
          i == 0 ? current_position_of_A[1] : current_position_of_B[1],
          i == 0 ? current_position_of_A[2] : current_position_of_B[2]}};

      const Block<3>& block =
          domain.blocks()[block_logical_coords[i].value().id.get_index()];
      const auto& grid_to_inertial_map =
          block.moving_mesh_grid_to_inertial_map();

      const auto inertial_position =
          grid_to_inertial_map(current_position, time, functions_of_time);

      for (size_t j = 0; j < 3; j++) {
        control_error[i * 3 + j] =
            inertial_position[0] - fot_inertial_positions[i * 3 + j];
      }
    }

    return control_error;
  }
};
}  // namespace control_system::ControlErrors
