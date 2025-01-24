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
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
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

  struct FrameToControl {
    using type = std::string;
    static constexpr Options::String help = {"Either 'Inertial' or 'Grid'."};
  };

  using options = tmpl::list<FrameToControl>;
  static constexpr Options::String help{
      "Computes the control error for the inertial centers of two objects. "
      "This should not take any options."};

  InertialCenters() = default;
  InertialCenters(std::string frame, const Options::Context& context = {})
      : frame_(std::move(frame)) {
    if (frame_ != "Inertial" and frame_ != "Grid") {
      PARSE_ERROR(context, "FrameToControl must be 'Inertial' or 'Grid'. Not "
                               << frame_);
    }
  }

  void pup(PUP::er& /*p*/) {}

  template <typename Metavariables, typename... TupleTags>
  DataVector operator()(const ::TimescaleTuner<true>& /*unused*/,
                        const Parallel::GlobalCache<Metavariables>& cache,
                        const double time,
                        const std::string& function_of_time_name,
                        const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
    const auto& function_of_time = functions_of_time.at(function_of_time_name);
    const DataVector fot_inertial_positions_dv =
        function_of_time->func(time)[0];

    using grid_center_A =
        control_system::QueueTags::Center<::domain::ObjectLabel::A,
                                          Frame::Grid>;
    using grid_center_B =
        control_system::QueueTags::Center<::domain::ObjectLabel::B,
                                          Frame::Grid>;
    using inertial_center_A =
        control_system::QueueTags::Center<::domain::ObjectLabel::A,
                                          Frame::Inertial>;
    using inertial_center_B =
        control_system::QueueTags::Center<::domain::ObjectLabel::B,
                                          Frame::Inertial>;

    const auto& measured_grid_position_of_A = get<grid_center_A>(measurements);
    const auto& measured_grid_position_of_B = get<grid_center_B>(measurements);
    const auto& measured_inertial_position_of_A =
        get<inertial_center_A>(measurements);
    const auto& measured_inertial_position_of_B =
        get<inertial_center_B>(measurements);

    tnsr::I<DataVector, 3, Frame::Grid> measured_grid_positions_tnsr{2_st, 0.0};
    for (size_t i = 0; i < 3; i++) {
      measured_grid_positions_tnsr.get(i)[0] = measured_grid_position_of_A[i];
      measured_grid_positions_tnsr.get(i)[1] = measured_grid_position_of_B[i];
    }

    tnsr::I<DataVector, 3, Frame::Inertial> fot_inertial_positions_tnsr{2_st,
                                                                        0.0};
    for (size_t i = 0; i < 3; i++) {
      fot_inertial_positions_tnsr.get(i)[0] = fot_inertial_positions_dv[i];
      fot_inertial_positions_tnsr.get(i)[1] = fot_inertial_positions_dv[3 + i];
    }

    const Domain<3>& domain = Parallel::get<domain::Tags::Domain<3>>(cache);

    // Get logical coords of measured centers
    const auto measured_logical_coords = block_logical_coordinates(
        domain, measured_grid_positions_tnsr, time, functions_of_time);
    ASSERT(alg::all_of(measured_logical_coords,
                       [](const auto& coord) { return coord.has_value(); }),
           "Measured centers are no longer in the domain. "
               << measured_grid_positions_tnsr);

    // Get logical coords of FoT centers
    const auto fot_logical_coords = block_logical_coordinates(
        domain, fot_inertial_positions_tnsr, time, functions_of_time);
    ASSERT(alg::all_of(fot_logical_coords,
                       [](const auto& coord) { return coord.has_value(); }),
           "FoT centers are no longer in the domain. "
               << fot_inertial_positions_tnsr);

    DataVector control_error{6, 0.0};

    // Compare the measured vs the FoT centers in either the grid or inertial
    // frame
    for (size_t i = 0; i < 2; i++) {
      if (frame_ == "Inertial") {
        // We already have the measured inertial positions, so just use them
        // directly to compare against the FoT inertial positions
        for (size_t j = 0; j < 3; j++) {
          control_error[i * 3 + j] =
              (i == 0 ? measured_inertial_position_of_A[j]
                      : measured_inertial_position_of_B[j]) -
              fot_inertial_positions_dv[i * 3 + j];
        }
      } else {
        // Need to conver the FoT inertial coords to the grid frame to compare
        // to the measured grid centers
        const auto& measured_grid_position =
            i == 0 ? measured_grid_position_of_A : measured_grid_position_of_B;
        const Block<3>& block =
            domain.blocks()[fot_logical_coords[i].value().id.get_index()];
        const auto& logical_to_grid_map =
            block.moving_mesh_logical_to_grid_map();

        const tnsr::I<double, 3, Frame::Grid> fot_grid_position =
            logical_to_grid_map(fot_logical_coords[i].value().data);

        for (size_t j = 0; j < 3; j++) {
          control_error[i * 3 + j] =
              measured_grid_position[j] - fot_grid_position.get(j);
        }
      }
    }

    return control_error;
  }

 private:
  std::string frame_{};
};
}  // namespace control_system::ControlErrors
