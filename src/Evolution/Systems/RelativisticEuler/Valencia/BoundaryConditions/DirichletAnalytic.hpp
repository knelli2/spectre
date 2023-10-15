// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <type_traits>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Evolution/Systems/RelativisticEuler/Valencia/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/Systems/RelativisticEuler/Valencia/ConservativeFromPrimitive.hpp"
#include "Evolution/Systems/RelativisticEuler/Valencia/Fluxes.hpp"
#include "Evolution/Systems/RelativisticEuler/Valencia/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct Time;
}  // namespace Tags
namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace domain::Tags
/// \endcond

namespace RelativisticEuler::Valencia::BoundaryConditions {
/*!
 * \brief Sets Dirichlet boundary conditions using the analytic solution or
 * analytic data.
 */
template <size_t Dim>
class DirichletAnalytic final : public BoundaryCondition<Dim> {
 public:
  using options = tmpl::list<>;
  inline const static std::string help{
      "DirichletAnalytic boundary conditions using either analytic solution or "
      "analytic data."};

  DirichletAnalytic() = default;
  DirichletAnalytic(DirichletAnalytic&&) = default;
  DirichletAnalytic& operator=(DirichletAnalytic&&) = default;
  DirichletAnalytic(const DirichletAnalytic&) = default;
  DirichletAnalytic& operator=(const DirichletAnalytic&) = default;
  ~DirichletAnalytic() override = default;

  explicit DirichletAnalytic(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, DirichletAnalytic);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  static constexpr evolution::BoundaryConditions::Type bc_type =
      evolution::BoundaryConditions::Type::Ghost;

  void pup(PUP::er& p) override;

  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags =
      tmpl::list<domain::Tags::Coordinates<Dim, Frame::Inertial>>;
  using dg_interior_primitive_variables_tags = tmpl::list<>;
  using dg_gridless_tags =
      tmpl::list<::Tags::Time, ::Tags::AnalyticSolutionOrData>;

  template <typename AnalyticSolutionOrData>
  std::optional<std::string> dg_ghost(
      const gsl::not_null<Scalar<DataVector>*> tilde_d,
      const gsl::not_null<Scalar<DataVector>*> tilde_tau,
      const gsl::not_null<tnsr::i<DataVector, Dim, Frame::Inertial>*> tilde_s,

      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          flux_tilde_d,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          flux_tilde_tau,
      const gsl::not_null<tnsr::Ij<DataVector, Dim, Frame::Inertial>*>
          flux_tilde_s,

      const gsl::not_null<Scalar<DataVector>*> lapse,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> shift,
      const gsl::not_null<tnsr::ii<DataVector, Dim, Frame::Inertial>*>
          spatial_metric,
      const gsl::not_null<Scalar<DataVector>*> rest_mass_density,
      const gsl::not_null<Scalar<DataVector>*> specific_internal_energy,
      const gsl::not_null<Scalar<DataVector>*> specific_enthalpy,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          spatial_velocity,
      const gsl::not_null<tnsr::II<DataVector, Dim, Frame::Inertial>*>
          inv_spatial_metric,

      const std::optional<
          tnsr::I<DataVector, Dim, Frame::Inertial>>& /*face_mesh_velocity*/,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& /*normal_covector*/,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& /*normal_vector*/,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& coords,
      const double time,
      const AnalyticSolutionOrData& analytic_solution_or_data) const {
    auto boundary_values = [&analytic_solution_or_data, &coords, &time]() {
      if constexpr (is_analytic_solution_v<AnalyticSolutionOrData>) {
        return analytic_solution_or_data.variables(
            coords, time,
            tmpl::list<hydro::Tags::RestMassDensity<DataVector>,
                       hydro::Tags::SpecificInternalEnergy<DataVector>,
                       hydro::Tags::SpecificEnthalpy<DataVector>,
                       hydro::Tags::Pressure<DataVector>,
                       hydro::Tags::SpatialVelocity<DataVector, Dim>,
                       hydro::Tags::LorentzFactor<DataVector>,
                       gr::Tags::SqrtDetSpatialMetric<DataVector>,
                       gr::Tags::Lapse<DataVector>,
                       gr::Tags::Shift<DataVector, Dim>,
                       gr::Tags::SpatialMetric<DataVector, Dim>,
                       gr::Tags::InverseSpatialMetric<DataVector, Dim>>{});
      } else {
        (void)time;
        return analytic_solution_or_data.variables(
            coords,
            tmpl::list<hydro::Tags::RestMassDensity<DataVector>,
                       hydro::Tags::SpecificInternalEnergy<DataVector>,
                       hydro::Tags::SpecificEnthalpy<DataVector>,
                       hydro::Tags::Pressure<DataVector>,
                       hydro::Tags::SpatialVelocity<DataVector, Dim>,
                       hydro::Tags::LorentzFactor<DataVector>,
                       gr::Tags::SqrtDetSpatialMetric<DataVector>,
                       gr::Tags::Lapse<DataVector>,
                       gr::Tags::Shift<DataVector, Dim>,
                       gr::Tags::SpatialMetric<DataVector, Dim>,
                       gr::Tags::InverseSpatialMetric<DataVector, Dim>>{});
      }
    }();

    *lapse = get<gr::Tags::Lapse<DataVector>>(boundary_values);
    *shift = get<gr::Tags::Shift<DataVector, Dim>>(boundary_values);
    *spatial_metric =
        get<gr::Tags::SpatialMetric<DataVector, Dim>>(boundary_values);
    *inv_spatial_metric =
        get<gr::Tags::InverseSpatialMetric<DataVector, Dim>>(boundary_values);
    *rest_mass_density =
        get<hydro::Tags::RestMassDensity<DataVector>>(boundary_values);
    *specific_internal_energy =
        get<hydro::Tags::SpecificInternalEnergy<DataVector>>(boundary_values);
    *specific_enthalpy =
        get<hydro::Tags::SpecificEnthalpy<DataVector>>(boundary_values);
    *spatial_velocity =
        get<hydro::Tags::SpatialVelocity<DataVector, Dim>>(boundary_values);

    const auto& pressure =
        get<hydro::Tags::Pressure<DataVector>>(boundary_values);
    const auto& lorentz_factor =
        get<hydro::Tags::LorentzFactor<DataVector>>(boundary_values);
    const auto& sqrt_det_spatial_metric =
        get<gr::Tags::SqrtDetSpatialMetric<DataVector>>(boundary_values);

    ConservativeFromPrimitive<Dim>::apply(
        tilde_d, tilde_tau, tilde_s, *rest_mass_density,
        *specific_internal_energy, *specific_enthalpy, pressure,
        *spatial_velocity, lorentz_factor, sqrt_det_spatial_metric,
        *spatial_metric);
    ComputeFluxes<Dim>::apply(flux_tilde_d, flux_tilde_tau, flux_tilde_s,
                              *tilde_d, *tilde_tau, *tilde_s, *lapse, *shift,
                              sqrt_det_spatial_metric, pressure,
                              *spatial_velocity);

    return {};
  }
};
}  // namespace RelativisticEuler::Valencia::BoundaryConditions
