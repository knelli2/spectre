// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/FixedHashMap.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/MaxNumberOfNeighbors.hpp"
#include "Domain/Structure/OrientationMapHelpers.hpp"
#include "Domain/Structure/TrimMap.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/ActiveGrid.hpp"
#include "Evolution/DgSubcell/GhostData.hpp"
#include "Evolution/DgSubcell/NeighborRdmpAndVolumeData.hpp"
#include "Evolution/DgSubcell/Projection.hpp"
#include "Evolution/DgSubcell/RdmpTci.hpp"
#include "Evolution/DgSubcell/RdmpTciData.hpp"
#include "Evolution/DgSubcell/SliceData.hpp"
#include "Evolution/DgSubcell/Tags/ActiveGrid.hpp"
#include "Evolution/DgSubcell/Tags/CellCenteredFlux.hpp"
#include "Evolution/DgSubcell/Tags/DataForRdmpTci.hpp"
#include "Evolution/DgSubcell/Tags/GhostDataForReconstruction.hpp"
#include "Evolution/DgSubcell/Tags/Mesh.hpp"
#include "Evolution/DgSubcell/Tags/TciStatus.hpp"
#include "Evolution/DiscontinuousGalerkin/InboxTags.hpp"
#include "Evolution/DiscontinuousGalerkin/Messages/BoundaryMessage.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarData.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/DiscontinuousGalerkin/Tags/NeighborMesh.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

namespace evolution::dg::subcell::Actions {
/*!
 * \brief Sets the local data from the relaxed discrete maximum principle
 * troubled-cell indicator and sends ghost zone data to neighboring elements.
 *
 * The action proceeds as follows:
 *
 * 1. Determine in which directions we have neighbors
 * 2. Slice the variables provided by GhostDataMutator to send to our neighbors
 *    for ghost zones
 * 3. Send the ghost zone data, appending the max/min for the TCI at the end of
 *    the `DataVector` we are sending.
 *
 * \warning This assumes the RDMP TCI data in the DataBox has been set, it does
 * not calculate it automatically. The reason is this way we can only calculate
 * the RDMP data when it's needed since computing it can be pretty expensive.
 *
 * Some notes:
 * - In the future we will need to send the cell-centered fluxes to do
 *   high-order FD without additional reconstruction being necessary.
 *
 * GlobalCache:
 * - Uses:
 *   - `ParallelComponent` proxy
 *
 * DataBox:
 * - Uses:
 *   - `domain::Tags::Mesh<Dim>`
 *   - `subcell::Tags::Mesh<Dim>`
 *   - `domain::Tags::Element<Dim>`
 *   - `Tags::TimeStepId`
 *   - `Tags::Next<Tags::TimeStepId>`
 *   - `subcell::Tags::ActiveGrid`
 *   - `System::variables_tag`
 *   - `subcell::Tags::DataForRdmpTci`
 * - Adds: nothing
 * - Removes: nothing
 * - Modifies:
 *   - `subcell::Tags::GhostDataForReconstruction<Dim>`
 */
template <size_t Dim, typename GhostDataMutator, bool LocalTimeStepping>
struct SendDataForReconstruction {
  using inbox_tags = tmpl::list<evolution::dg::Tags::BoundaryMessageInbox<Dim>>;

  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent,
            typename Metavariables>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    static_assert(
        not LocalTimeStepping,
        "DG-subcell does not yet support local time stepping. The "
        "reconstruction data must be sent using dense output sometimes, and "
        "not at all other times. However, the data for the RDMP TCI should be "
        "sent along with the data for reconstruction each time.");

    ASSERT(db::get<Tags::ActiveGrid>(box) == ActiveGrid::Subcell,
           "The SendDataForReconstruction action can only be called when "
           "Subcell is the active scheme.");
    using flux_variables = typename Metavariables::system::flux_variables;

    // QUESTION: Do we need this?
    // db::mutate<Tags::GhostDataForReconstruction<Dim>>(
    //     make_not_null(&box), [](const auto ghost_data_ptr) {
    //       // Clear the previous neighbor data and add current local data
    //       ghost_data_ptr->clear();
    //     });

    const Mesh<Dim>& dg_mesh = db::get<::domain::Tags::Mesh<Dim>>(box);
    const Mesh<Dim>& subcell_mesh = db::get<Tags::Mesh<Dim>>(box);
    const Element<Dim>& element = db::get<::domain::Tags::Element<Dim>>(box);
    const size_t ghost_zone_size =
        Metavariables::SubcellOptions::ghost_zone_size(box);
    DirectionMap<Dim, bool> directions_to_slice{};
    for (const auto& direction_neighbors : element.neighbors()) {
      if (direction_neighbors.second.size() == 0) {
        directions_to_slice[direction_neighbors.first] = false;
      } else {
        directions_to_slice[direction_neighbors.first] = true;
      }
    }

    // Optimization note: could save a copy+allocation if we moved
    // all_sliced_data when possible before sending.
    //
    // Note: RDMP size doesn't help here since we need to slice data after
    // anyway, so no way to save an allocation through that.
    const auto& cell_centered_flux =
        db::get<Tags::CellCenteredFlux<flux_variables, Dim>>(box);
    DataVector volume_data_to_slice = db::mutate_apply(
        GhostDataMutator{}, make_not_null(&box),
        cell_centered_flux.has_value() ? cell_centered_flux.value().size()
                                       : 0_st);
    if (cell_centered_flux.has_value()) {
      std::copy(
          cell_centered_flux.value().data(),
          std::next(
              cell_centered_flux.value().data(),
              static_cast<std::ptrdiff_t>(cell_centered_flux.value().size())),
          std::next(
              volume_data_to_slice.data(),
              static_cast<std::ptrdiff_t>(volume_data_to_slice.size() -
                                          cell_centered_flux.value().size())));
    }

    auto& subcell_neighbor_data_pair = db::get_mutable_reference<
        evolution::dg::subcell::Tags::NeighborData<Dim>>(make_not_null(&box));
    const size_t neighbor_data_index = subcell_neighbor_data_pair.first;
    DirectionMap<Dim, DataVector>& all_sliced_data =
        subcell_neighbor_data_pair.second[neighbor_data_index].value();

    auto& receiver_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const RdmpTciData& rdmp_tci_data = db::get<Tags::DataForRdmpTci>(box);
    const TimeStepId& time_step_id = db::get<::Tags::TimeStepId>(box);
    const TimeStepId& next_time_step_id = [&box]() {
      if (LocalTimeStepping) {
        return db::get<::Tags::Next<::Tags::TimeStepId>>(box);
      } else {
        return db::get<::Tags::TimeStepId>(box);
      }
    }();

    const size_t rdmp_size = rdmp_tci_data.max_variables_values.size() +
                             rdmp_tci_data.min_variables_values.size();

    all_sliced_data =
        slice_data(volume_data_to_slice, subcell_mesh.extents(),
                   ghost_zone_size, directions_to_slice, rdmp_size);

    const int tci_decision =
        db::get<evolution::dg::subcell::Tags::TciDecision>(box);
    // Compute and send actual variables
    for (const auto& [direction, neighbors_in_direction] :
         element.neighbors()) {
      const auto& orientation = neighbors_in_direction.orientation();
      const auto direction_from_neighbor = orientation(direction.opposite());
      ASSERT(neighbors_in_direction.size() == 1,
             "AMR is not yet supported when using DG-subcell. Note that this "
             "condition could be relaxed to support AMR only where the "
             "evolution is using DG without any changes to subcell.");

      for (const ElementId<Dim>& neighbor : neighbors_in_direction) {
        const DataVector& sliced_data_in_direction =
            all_sliced_data.at(direction);

        using BoundaryMessage = evolution::dg::BoundaryMessage<Dim>;

        BoundaryMessage* message = nullptr;

        const bool aligned_with_neighbor = orientation.is_aligned();
        const size_t subcell_size = sliced_data_in_direction.size() + rdmp_size;
        constexpr size_t dg_size = 0;

        // If we are aligned with neighbor, we make a non-owning message. But if
        // we need to align our data for our neighbor, we create an owning
        // message to send because we must preserve our own local_mortar_data
        if (LIKELY(aligned_with_neighbor)) {
          message = new BoundaryMessage{};
          message->owning = false;
        } else {
          const size_t total_bytes_with_data =
              BoundaryMessage::total_bytes_with_data(subcell_size, dg_size);

          // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
          message = reinterpret_cast<BoundaryMessage*>(
              BoundaryMessage::operator new(total_bytes_with_data));
          message->owning = true;
        }

        // Set everything except the FD pointer
        message->enable_if_disabled = false;
        message->sender_node = Parallel::my_node<size_t>(cache);
        message->sender_core = Parallel::my_proc<size_t>(cache);
        message->tci_status = tci_decision;
        message->current_time_step_id = time_step_id;
        message->next_time_step_id = next_time_step_id;
        message->neighbor_direction = direction_from_neighbor;
        message->element_id = element.id();
        message->volume_or_ghost_mesh = subcell_mesh;
        message->interface_mesh = dg_mesh.slice_away(direction.dimension());
        message->subcell_ghost_data_size = subcell_size;
        message->dg_flux_data_size = dg_size;
        // We don't send DG data here
        message->dg_flux_data = nullptr;

        // If we are aligned, then the subcell pointer is to the data stored in
        // the Tags::Neighbor data tag
        if (LIKELY(aligned_with_neighbor)) {
          message->subcell_ghost_data =
              const_cast<double*>(sliced_data_in_direction.data());
        } else {
          // If we need to orient, we avoid extra allocations/copies by
          // orienting directly into the message buffer we allocated above.

          // double* + 1 == char* + 8 because double* is 8 bytes
          // Place subcell data right after dg pointer
          message->subcell_ghost_data =
              // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
              reinterpret_cast<double*>(std::addressof(message->dg_flux_data))
              // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
              + 1;

          std::array<size_t, Dim> slice_extents{};
          for (size_t d = 0; d < Dim; ++d) {
            gsl::at(slice_extents, d) = subcell_mesh.extents(d);
          }
          gsl::at(slice_extents, direction.dimension()) = ghost_zone_size;

          // Need a view so we only get the subcell data and not the rdmp data
          DataVector subcell_data_to_send_view{
              message->subcell_ghost_data,
              message->subcell_ghost_data_size - rdmp_size};

          orient_variables(make_not_null(&subcell_data_to_send_view),
                           sliced_data_in_direction, Index<Dim>{slice_extents},
                           orientation);
        }

        // Copy rdmp data to end of message->subcell_ghost_data. The pointer
        // always points to the proper place we want. If we are oriented, then
        // the pointer is to the sliced_data in the Tags::NeighborData tag. If
        // we aren't oriented, then the pointer is to the buffer we just
        // allocated.
        DataVector subcell_data_to_send{message->subcell_ghost_data,
                                        message->subcell_ghost_data_size};
        std::copy(
            rdmp_tci_data.max_variables_values.cbegin(),
            rdmp_tci_data.max_variables_values.cend(),
            std::prev(subcell_data_to_send.end(), static_cast<int>(rdmp_size)));
        std::copy(rdmp_tci_data.min_variables_values.cbegin(),
                  rdmp_tci_data.min_variables_values.cend(),
                  std::prev(subcell_data_to_send.end(),
                            static_cast<int>(
                                rdmp_tci_data.min_variables_values.size())));

        Parallel::receive_data<evolution::dg::Tags::BoundaryMessageInbox<Dim>>(
            receiver_proxy[neighbor], message);
      }
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

/*!
 * \brief Receive the subcell data from our neighbor, and accumulate the data
 * from the relaxed discrete maximum principle troubled-cell indicator.
 *
 * Note:
 * - Since we only care about the min/max over all neighbors and ourself at the
 *   past time, we accumulate all data immediately into the `RdmpTciData`.
 * - If the neighbor is using DG and therefore sends boundary correction data
 *   then that is added into the `evolution::dg::Tags::MortarData` tag
 * - The next `TimeStepId` is recorded, but we do not yet support local time
 *   stepping.
 * - This action will never care about what variables are sent for
 *   reconstruction. It is only responsible for receiving the data and storing
 *   it in the `NeighborData`.
 *
 * GlobalCache:
 * -Uses: nothing
 *
 * DataBox:
 * - Uses:
 *   - `domain::Tags::Element<Dim>`
 *   - `Tags::TimeStepId`
 *   - `domain::Tags::Mesh<Dim>`
 *   - `subcell::Tags::Mesh<Dim>`
 *   - `domain::Tags::Element<Dim>`
 *   - `Tags::Next<Tags::TimeStepId>`
 *   - `subcell::Tags::ActiveGrid`
 *   - `System::variables_tag`
 * - Adds: nothing
 * - Removes: nothing
 * - Modifies:
 *   - `subcell::Tags::GhostDataForReconstruction<Dim>`
 *   - `subcell::Tags::DataForRdmpTci`
 *   - `evolution::dg::Tags::MortarData`
 *   - `evolution::dg::Tags::MortarNextTemporalId`
 */
template <size_t Dim>
struct ReceiveDataForReconstruction {
  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent,
            typename Metavariables>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const Element<Dim>& element = db::get<::domain::Tags::Element<Dim>>(box);
    const auto number_of_expected_messages = element.neighbors().size();
    if (UNLIKELY(number_of_expected_messages == 0)) {
      // We have no neighbors, so just continue without doing any work
      return {Parallel::AlgorithmExecution::Continue, std::nullopt};
    }

    using ::operator<<;
    using Key = std::pair<Direction<Dim>, ElementId<Dim>>;
    using BoundaryMessage = evolution::dg::BoundaryMessage<Dim>;
    using MessageType = std::unique_ptr<BoundaryMessage>;
    using MapType = FixedHashMap<maximum_number_of_neighbors(Dim), Key,
                                 MessageType, boost::hash<Key>>;
    using InboxType = std::map<TimeStepId, MapType>;

    const auto& current_time_step_id = db::get<::Tags::TimeStepId>(box);
    InboxType& inbox = tuples::get<
        evolution::dg::Tags::BoundaryMessageInbox<Metavariables::volume_dim>>(
        inboxes);
    const auto received = inbox.find(current_time_step_id);
    // Check we have at least some data from correct time, and then check that
    // we have received all data
    if (received == inbox.end() or
        received->second.size() != number_of_expected_messages) {
      return {Parallel::AlgorithmExecution::Retry, std::nullopt};
    }

    MapType& received_data = received->second;

    const Mesh<Dim>& subcell_mesh = db::get<Tags::Mesh<Dim>>(box);

    db::mutate<Tags::GhostDataForReconstruction<Dim>, Tags::DataForRdmpTci,
               evolution::dg::Tags::MortarData<Dim>,
               evolution::dg::Tags::MortarNextTemporalId<Dim>,
               evolution::dg::Tags::NeighborMesh<Dim>,
               evolution::dg::Tags::BoundaryMessageFromInbox<Dim>,
               evolution::dg::subcell::Tags::NeighborTciDecisions<Dim>>(
        make_not_null(&box),
        [&current_time_step_id, &element,
         ghost_zone_size = Metavariables::SubcellOptions::ghost_zone_size(box),
         &received_data, &subcell_mesh](
            const gsl::not_null<FixedHashMap<
                maximum_number_of_neighbors(Dim),
                std::pair<Direction<Dim>, ElementId<Dim>>, GhostData,
                boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>*>
                ghost_data_ptr,
            const gsl::not_null<RdmpTciData*> rdmp_tci_data_ptr,
            const gsl::not_null<std::unordered_map<
                Key, evolution::dg::MortarData<Dim>, boost::hash<Key>>*>
                mortar_data,
            const gsl::not_null<
                std::unordered_map<Key, TimeStepId, boost::hash<Key>>*>
                mortar_next_time_step_id,
            const gsl::not_null<FixedHashMap<
                maximum_number_of_neighbors(Dim),
                std::pair<Direction<Dim>, ElementId<Dim>>, Mesh<Dim>,
                boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>*>
                neighbor_mesh,
            const gsl::not_null<
                std::unordered_map<Key, MessageType, boost::hash<Key>>*>
                boundary_messages_in_data_box,
            const auto neighbor_tci_decisions) {
          // Remove neighbor meshes for neighbors that don't exist anymore
          domain::remove_nonexistent_neighbors(neighbor_mesh, element);

          // Get the next time step id, and also the fluxes data if the neighbor
          // is doing DG.
          for (std::pair<const Key, MessageType>& received_mortar_data :
               received_data) {
            const Key& mortar_id = received_mortar_data.first;
            MessageType& boundary_message_ptr = received_mortar_data.second;

            try {
              mortar_next_time_step_id->at(mortar_id) =
                  boundary_message_ptr->next_time_step_id;
            } catch (std::exception& e) {
              ERROR("Failed retrieving the MortarId: ("
                    << mortar_id.first << ',' << mortar_id.second
                    << ") from the mortar_next_time_step_id. Got exception: "
                    << e.what());
            }

            if (boundary_message_ptr->dg_flux_data_size != 0) {
              ASSERT(boundary_message_ptr->dg_flux_data != nullptr,
                     "In ReceiveDataForReconstruction, boundary message at "
                     "MortarId: ("
                         << mortar_id.first << "," << mortar_id.second
                         << ") has non-zero size for dg flux data, but has a "
                            "nullptr for the actual data.");

              auto& neighbor_mortar_data =
                  mortar_data->at(mortar_id).neighbor_mortar_data();

              if (UNLIKELY(not neighbor_mortar_data.has_value())) {
                neighbor_mortar_data = std::pair<Mesh<Dim - 1>, DataVector>{};
              }

              neighbor_mortar_data.value().first =
                  boundary_message_ptr->interface_mesh;
              neighbor_mortar_data.value().second.set_data_ref(
                  boundary_message_ptr->dg_flux_data,
                  boundary_message_ptr->dg_flux_data_size);
            }
            // Set new neighbor mesh
            neighbor_mesh->insert_or_assign(
                mortar_id, boundary_message_ptr->volume_or_ghost_mesh);
          }

          // ASSERT(ghost_data_ptr->empty(),
          //        "Should have no elements in the neighbor data when "
          //        "receiving neighbor data");
          const size_t number_of_rdmp_vars =
              rdmp_tci_data_ptr->max_variables_values.size();
          ASSERT(rdmp_tci_data_ptr->min_variables_values.size() ==
                     number_of_rdmp_vars,
                 "The number of RDMP variables for which we have a maximum "
                 "and minimum should be the same, but we have "
                     << number_of_rdmp_vars << " for the max and "
                     << rdmp_tci_data_ptr->min_variables_values.size()
                     << " for the min.");

          for (const auto& [direction, neighbors_in_direction] :
               element.neighbors()) {
            for (const auto& neighbor : neighbors_in_direction) {
              Key mortar_id{direction, neighbor};
              MessageType& boundary_message_ptr = received_data[mortar_id];
              // ASSERT(neighbor_data_ptr->count(mortar_id) == 0,
              //        "Found neighbor already inserted in direction "
              //            << direction << " with ElementId " << neighbor);
              ASSERT(boundary_message_ptr->subcell_ghost_data_size != 0 and
                         boundary_message_ptr->subcell_ghost_data != nullptr,
                     "Received subcell data message that does not contain any "
                     "actual subcell data for reconstruction.");
              // Collect the max/min of u(t^n) for the RDMP as we receive data.
              // This reduces the memory footprint.

              DataVector received_subcell_data{
                  boundary_message_ptr->subcell_ghost_data,
                  boundary_message_ptr->subcell_ghost_data_size};

              evolution::dg::subcell::insert_neighbor_rdmp_and_volume_data(
                  rdmp_tci_data_ptr, ghost_data_ptr, received_subcell_data,
                  number_of_rdmp_vars, mortar_id, neighbor_mesh->at(mortar_id),
                  element, subcell_mesh, ghost_zone_size);
              ASSERT(neighbor_tci_decisions->contains(mortar_id),
                     "The NeighorTciDecisions should contain the neighbor ("
                         << mortar_id.first << ", " << mortar_id.second
                         << ") but doesn't");
              neighbor_tci_decisions->at(mortar_id) =
                  boundary_message_ptr->tci_status;

              // Now we move the boundary message from the inbox to the DataBox.
              // We do this because if the mortar data came from another node,
              // then the unique_ptr actually owns the mortar data. It isn't
              // owned by another DataBox. So if we only copy the fd/dg pointer
              // here, but then delete the BoundaryMessage pointer from the
              // inbox, then the mortar data owned by that unique_ptr will also
              // be deleted. These pointers don't need to be removed from the
              // DataBox because once they get overridden here, the previous
              // memory gets freed
              boundary_messages_in_data_box->at(mortar_id) =
                  std::move(boundary_message_ptr);
            }
          }
        });

    // Now we can erase the inbox
    inbox.erase(current_time_step_id);
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace evolution::dg::subcell::Actions
