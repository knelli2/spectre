// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <boost/functional/hash.hpp>
#include <cstddef>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/FixedHashMap.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/MaxNumberOfNeighbors.hpp"
#include "Evolution/DgSubcell/GhostData.hpp"
#include "Evolution/DgSubcell/NeighborReconstructedFaceSolution.hpp"
#include "Evolution/DgSubcell/RdmpTciData.hpp"
#include "Evolution/DgSubcell/Tags/DataForRdmpTci.hpp"
#include "Evolution/DgSubcell/Tags/GhostDataForReconstruction.hpp"
#include "Evolution/DiscontinuousGalerkin/Messages/BoundaryMessage.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct VolumeDouble : db::SimpleTag {
  using type = double;
};

template <size_t Dim>
using GhostDataMap =
    FixedHashMap<maximum_number_of_neighbors(Dim),
                 std::pair<Direction<Dim>, ElementId<Dim>>,
                 evolution::dg::subcell::GhostData,
                 boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>;
template <size_t Dim>
using NeighborReconstructionMap =
    FixedHashMap<maximum_number_of_neighbors(Dim),
                 std::pair<Direction<Dim>, ElementId<Dim>>, DataVector,
                 boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>;

template <size_t Dim>
using BoundaryMessage = evolution::dg::BoundaryMessage<Dim>;

template <size_t Dim>
using MortarData = std::unique_ptr<BoundaryMessage<Dim>>;

template <size_t Dim>
using MortarDataMap =
    FixedHashMap<maximum_number_of_neighbors(Dim),
                 std::pair<Direction<Dim>, ElementId<Dim>>, MortarData<Dim>,
                 boost::hash<std::pair<Direction<Dim>, ElementId<Dim>>>>;

template <size_t Dim>
struct Metavariables {
  static constexpr size_t volume_dim = Dim;
  struct SubcellOptions {
    struct DgComputeSubcellNeighborPackagedData {
      template <typename DbTagsList>
      static NeighborReconstructionMap<Dim> apply(
          const db::DataBox<DbTagsList>& box,
          const std::vector<
              std::pair<Direction<volume_dim>, ElementId<volume_dim>>>&
              mortars_to_reconstruct_to) {
        const GhostDataMap<Dim>& ghost_data = db::get<
            evolution::dg::subcell::Tags::GhostDataForReconstruction<Dim>>(box);

        // We just simply copy over the data sent since it doesn't actually
        // matter what we fill the packaged data with in the test, just that
        // this function is called and that we can retrieve the correct data
        // from the stored NeighborData.
        NeighborReconstructionMap<Dim> neighbor_package_data{};
        for (const auto& mortar_id : mortars_to_reconstruct_to) {
          neighbor_package_data[mortar_id] =
              ghost_data.at(mortar_id).neighbor_ghost_data_for_reconstruction();
        }
        return neighbor_package_data;
      }
    };
  };
};

template <size_t Dim>
BoundaryMessage<Dim>* create_message(
    const Mesh<Dim>& volume_mesh, const Mesh<Dim - 1>& face_mesh,
    const DataVector& subcell_data,
    const std::optional<DataVector>& dg_flux_data, const int tci_status) {
  BoundaryMessage<Dim>* message = new BoundaryMessage<Dim>();
  message->volume_or_ghost_mesh = volume_mesh;
  message->interface_mesh = face_mesh;
  message->tci_status = tci_status;

  message->subcell_ghost_data_size = subcell_data.size();
  message->subcell_ghost_data = const_cast<double*>(subcell_data.data());

  message->dg_flux_data_size = dg_flux_data.value_or(DataVector{}).size();
  message->dg_flux_data = dg_flux_data.has_value()
                              ? const_cast<double*>(dg_flux_data.value().data())
                              : nullptr;

  return message;
}

template <size_t Dim>
void test() {
  CAPTURE(Dim);
  using metavars = Metavariables<Dim>;
  const std::pair self_id{Direction<Dim>::lower_xi(),
                          ElementId<Dim>::external_boundary_id()};

  evolution::dg::subcell::RdmpTciData rdmp_tci_data{};
  rdmp_tci_data.max_variables_values = DataVector{1.0, 2.0};
  rdmp_tci_data.min_variables_values = DataVector{-2.0, 0.1};
  GhostDataMap<Dim> neighbor_data_map{};
  auto box = db::create<
      tmpl::list<evolution::dg::subcell::Tags::GhostDataForReconstruction<Dim>,
                 evolution::dg::subcell::Tags::DataForRdmpTci, VolumeDouble>>(
      std::move(neighbor_data_map), std::move(rdmp_tci_data), 2.5);

  std::pair<const TimeStepId, MortarDataMap<Dim>> mortar_data_from_neighbors{};
  std::array<DataVector, Dim> fd_recons_and_rdmp_data_array{};
  std::array<DataVector, Dim> dg_recons_and_rdmp_data_array{};
  for (size_t d = 0; d < Dim; ++d) {
    const Mesh<Dim> dg_volume_mesh{2 + 2 * Dim, Spectral::Basis::Legendre,
                                   Spectral::Quadrature::GaussLobatto};
    const Mesh<Dim> fd_volume_mesh{2 + 2 * Dim + 1,
                                   Spectral::Basis::FiniteDifference,
                                   Spectral::Quadrature::CellCentered};
    const Mesh<Dim - 1> dg_face_mesh{2 + 2 * Dim, Spectral::Basis::Legendre,
                                     Spectral::Quadrature::GaussLobatto};
    const Mesh<Dim - 1> fd_face_mesh{2 + 2 * Dim + 1,
                                     Spectral::Basis::FiniteDifference,
                                     Spectral::Quadrature::CellCentered};
    DataVector& fd_recons_and_rdmp_data =
        gsl::at(fd_recons_and_rdmp_data_array, d);
    DataVector& dg_recons_and_rdmp_data =
        gsl::at(dg_recons_and_rdmp_data_array, d);
    fd_recons_and_rdmp_data = DataVector{2 * Dim + 1 + 4, 4.0};
    dg_recons_and_rdmp_data = DataVector{2 * Dim + 1 + 4, 7.0};
    for (size_t i = 0; i < 4; ++i) {
      dg_recons_and_rdmp_data[2 * Dim + 1 + i] =
          (i > 1 ? -1.0 : 1.0) * (d + 1.0) * 7.0 * (i + 5.0);
      fd_recons_and_rdmp_data[2 * Dim + 1 + i] =
          (i > 1 ? -1.0 : 1.0) * (d + 1.0) * 7.0 * (i + 50.0);
    }
    DataVector dg_flux_data(2 * Dim + 1);
    if (d % 2 == 0) {
      mortar_data_from_neighbors.second[std::pair{
          Direction<Dim>{d, Side::Upper}, ElementId<Dim>{2 * d}}] =
          std::unique_ptr<BoundaryMessage<Dim>>(create_message(
              dg_volume_mesh, dg_face_mesh, dg_recons_and_rdmp_data,
              std::optional{dg_flux_data}, 1));
      mortar_data_from_neighbors.second[std::pair{
          Direction<Dim>{d, Side::Lower}, ElementId<Dim>{2 * d + 1}}] =
          std::unique_ptr<BoundaryMessage<Dim>>(
              create_message(fd_volume_mesh, fd_face_mesh,
                             fd_recons_and_rdmp_data, std::nullopt, 2));
    } else {
      mortar_data_from_neighbors.second[std::pair{
          Direction<Dim>{d, Side::Lower}, ElementId<Dim>{2 * d}}] =
          std::unique_ptr<BoundaryMessage<Dim>>(create_message(
              dg_volume_mesh, dg_face_mesh, dg_recons_and_rdmp_data,
              std::optional{dg_flux_data}, 3));
      mortar_data_from_neighbors.second[std::pair{
          Direction<Dim>{d, Side::Upper}, ElementId<Dim>{2 * d + 1}}] =
          std::unique_ptr<BoundaryMessage<Dim>>(
              create_message(fd_volume_mesh, fd_face_mesh,
                             fd_recons_and_rdmp_data, std::nullopt, 4));
    }
  }

  evolution::dg::subcell::neighbor_reconstructed_face_solution<metavars>(
      make_not_null(&box), make_not_null(&mortar_data_from_neighbors));

  for (size_t d = 0; d < Dim; ++d) {
    CAPTURE(d);
    const bool d_is_odd = (d % 2 != 0);
    const std::pair id{Direction<Dim>{d, d_is_odd ? Side::Upper : Side::Lower},
                       ElementId<Dim>{2 * d + 1}};
    CAPTURE(id);

    const MortarDataMap<Dim>& neighbor_mortar_map =
        mortar_data_from_neighbors.second;
    REQUIRE(neighbor_mortar_map.contains(id));

    const MortarData<Dim>& boundary_message = neighbor_mortar_map.at(id);

    REQUIRE(boundary_message->subcell_ghost_data_size != 0);
    REQUIRE(boundary_message->subcell_ghost_data != nullptr);
    REQUIRE(boundary_message->dg_flux_data_size != 0);
    REQUIRE(boundary_message->dg_flux_data != nullptr);

    const DataVector dg_flux_data{boundary_message->dg_flux_data,
                                  boundary_message->dg_flux_data_size};
    const DataVector expected_dg_flux_data{
        boundary_message->subcell_ghost_data,
        boundary_message->subcell_ghost_data_size - 4};
    CHECK(dg_flux_data == expected_dg_flux_data);
    CHECK(boundary_message->tci_status == (d_is_odd ? 4 : 2));
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Subcell.NeighborReconstructedFaceSolution",
                  "[Evolution][Unit]") {
  test<1>();
  test<2>();
  test<3>();
}
