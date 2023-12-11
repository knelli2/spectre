// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/Protocols/FirstOrderSystem.hpp"
#include "Elliptic/Systems/Xcts/FluxesAndSources.hpp"
#include "Elliptic/Systems/Xcts/Geometry.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags/Conformal.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace Xcts {

/*!
 * \brief The Extended Conformal Thin Sandwich (XCTS) decomposition of the
 * Einstein constraint equations, formulated as a set of coupled first-order
 * partial differential equations
 *
 * See \ref Xcts for details on the XCTS equations.
 * The system can be formulated in terms of these fluxes and sources (see
 * `elliptic::protocols::FirstOrderSystem`):
 *
 * \f{align}
 * F^i_\psi ={} &\bar{\gamma}^{ij} \partial_j \psi \\
 * S_\psi ={} &-\bar{\Gamma}^i_{ij} F^j_\psi
 * + \frac{1}{8}\psi\bar{R} + \frac{1}{12}\psi^5 K^2
 * - \frac{1}{8}\psi^{-7}\bar{A}^2 - 2\pi\psi^5\rho
 * \f}
 *
 * for the Hamiltonian constraint,
 *
 * \f{align}
 * F^i_{\alpha\psi} ={} &\bar{\gamma}^{ij} \partial_j \alpha\psi \\
 * S_{\alpha\psi} ={} &-\bar{\Gamma}^i_{ij} F^j_{\alpha\psi}
 * + \alpha\psi \left(\frac{7}{8}\psi^{-8} \bar{A}^2
 * + \frac{5}{12} \psi^4 K^2 + \frac{1}{8}\bar{R}
 * + 2\pi\psi^4\left(\rho + 2S\right) \right) \\
 * &- \psi^5\partial_t K + \psi^5\left(\beta^i\bar{D}_i K
 * + \beta_\mathrm{background}^i\bar{D}_i K\right)
 * \f}
 *
 * for the lapse equation, and
 *
 * \f{align}
 * F^{ij}_\beta ={} &\left(\bar{L}\beta\right)^{ij} = \bar{\nabla}^i \beta^j
 * + \bar{\nabla}^j \beta^i
 * - \frac{2}{3} \bar{\gamma}^{ij} \bar{\nabla}_k \beta^k \\
 * S^i_\beta ={} &-\bar{\Gamma}^j_{jk} F^{ik}_\beta
 * - \bar{\Gamma}^i_{jk} F^{jk}_\beta
 * + \left(F^{ij}_\beta
 * + \left(\bar{L}\beta_\mathrm{background}\right)^{ij} - \bar{u}^{ij}\right)
 * \bar{\gamma}_{jk} \left(\frac{F^k_{\alpha\psi}}{\alpha\psi}
 * - 7 \frac{F^k_\psi}{\psi}\right) \\
 * &- \bar{D}_j\left(\left(\bar{L}\beta_\mathrm{background}\right)^{ij}
 * - \bar{u}^{ij}\right)
 * + \frac{4}{3}\frac{\alpha\psi}{\psi}\bar{D}^i K
 * + 16\pi\left(\alpha\psi\right)\psi^3 S^i
 * \f}
 *
 * for the momentum constraint, with
 *
 * \f{align}
 * \bar{A}^{ij} ={} &\frac{\psi^7}{2\alpha\psi}\left(
 * \left(\bar{L}\beta\right)^{ij} +
 * \left(\bar{L}\beta_\mathrm{background}\right)^{ij} - \bar{u}^{ij} \right)
 * \f}
 *
 * and all \f$f_\alpha=0\f$.
 *
 * Note that the symbol \f$\beta\f$ in the equations above means
 * \f$\beta_\mathrm{excess}\f$. The full shift is \f$\beta_\mathrm{excess} +
 * \beta_\mathrm{background}\f$. See `Xcts::Tags::ShiftBackground` and
 * `Xcts::Tags::ShiftExcess` for details on this split. Also note that the
 * background shift is degenerate with \f$\bar{u}\f$ so we treat the quantity
 * \f$\left(\bar{L}\beta_\mathrm{background}\right)^{ij} - \bar{u}^{ij}\f$ as a
 * single background field (see
 * `Xcts::Tags::LongitudinalShiftBackgroundMinusDtConformalMetric`). The
 * covariant divergence of this quantity w.r.t. the conformal metric is also a
 * background field.
 *
 * \par Solving a subset of equations:
 * This system allows you to select a subset of `Xcts::Equations` so you don't
 * have to solve for all variables if some are analytically known. Specify the
 * set of enabled equations as the first template parameter. The set of required
 * background fields depends on your choice of equations.
 *
 * \par Conformal background geometry:
 * The equations simplify significantly if the conformal metric is flat
 * ("conformal flatness") and in Cartesian coordinates. In this case you can
 * specify `Xcts::Geometry::FlatCartesian` as the second template parameter so
 * computations are optimized for a flat background geometry and you don't have
 * to supply geometric background fields. Else, specify
 * `Xcts::Geometry::Curved`.
 *
 * \par Conformal matter scale:
 * The matter source terms in the XCTS equations have the known defect that they
 * can spoil uniqueness of the solutions. See e.g. \cite Baumgarte2006ug for a
 * detailed study. To cure this defect one can conformally re-scale the matter
 * source terms as \f$\bar{\rho}=\psi^n\rho\f$, \f$\bar{S}=\psi^n S\f$ and
 * \f$\bar{S^i}=\psi^n S^i\f$ and treat the re-scaled fields as
 * freely-specifyable background data for the XCTS equations. You can select the
 * `ConformalMatterScale` \f$n\f$ as the third template parameter. Common
 * choices are \f$n=0\f$ for vacuum systems where the matter sources are
 * irrelevant, \f$n=6\f$ as suggested in \cite Foucart2008qt or \f$n=8\f$ as
 * suggested in \cite Baumgarte2006ug.
 */
template <Equations EnabledEquations, Geometry ConformalGeometry,
          int ConformalMatterScale>
struct FirstOrderSystem
    : tt::ConformsTo<elliptic::protocols::FirstOrderSystem> {
  static constexpr Equations enabled_equations = EnabledEquations;
  static constexpr Geometry conformal_geometry = ConformalGeometry;
  static constexpr int conformal_matter_scale = ConformalMatterScale;
  static constexpr size_t volume_dim = 3;

  using primal_fields = tmpl::flatten<tmpl::list<
      Tags::ConformalFactorMinusOne<DataVector>,
      tmpl::conditional_t<
          EnabledEquations == Equations::HamiltonianAndLapse or
              EnabledEquations == Equations::HamiltonianLapseAndShift,
          Tags::LapseTimesConformalFactorMinusOne<DataVector>, tmpl::list<>>,
      tmpl::conditional_t<
          EnabledEquations == Equations::HamiltonianLapseAndShift,
          Tags::ShiftExcess<DataVector, 3, Frame::Inertial>, tmpl::list<>>>>;

  // As fluxes we use the gradients with raised indices for the Hamiltonian and
  // lapse equation, and the longitudinal shift excess for the momentum
  // constraint. The gradient fluxes don't have symmetries and no particular
  // meaning so we use the standard `Flux` tags, but for the symmetric
  // longitudinal shift we use the corresponding symmetric tag.
  using primal_fluxes = tmpl::flatten<tmpl::list<
      ::Tags::Flux<Tags::ConformalFactorMinusOne<DataVector>, tmpl::size_t<3>,
                   Frame::Inertial>,
      tmpl::conditional_t<
          EnabledEquations == Equations::HamiltonianAndLapse or
              EnabledEquations == Equations::HamiltonianLapseAndShift,
          ::Tags::Flux<Tags::LapseTimesConformalFactorMinusOne<DataVector>,
                       tmpl::size_t<3>, Frame::Inertial>,
          tmpl::list<>>,
      tmpl::conditional_t<
          EnabledEquations == Equations::HamiltonianLapseAndShift,
          Tags::LongitudinalShiftExcess<DataVector, 3, Frame::Inertial>,
          tmpl::list<>>>>;

  using background_fields = tmpl::flatten<tmpl::list<
      // Quantities for Hamiltonian constraint
      gr::Tags::Conformal<gr::Tags::EnergyDensity<DataVector>,
                          ConformalMatterScale>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>,
      tmpl::conditional_t<
          ConformalGeometry == Geometry::Curved,
          tmpl::list<
              Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>,
              Tags::ConformalRicciScalar<DataVector>,
              Tags::ConformalChristoffelContracted<DataVector, 3,
                                                   Frame::Inertial>,
              ::Tags::deriv<
                  Tags::ConformalMetric<DataVector, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial>>,
          tmpl::list<>>,
      tmpl::conditional_t<
          EnabledEquations ==
              Equations::Hamiltonian,
          Tags::LongitudinalShiftMinusDtConformalMetricOverLapseSquare<
              DataVector>,
          tmpl::list<>>,
      // Additional quantities for lapse equation
      tmpl::conditional_t<
          EnabledEquations == Equations::HamiltonianAndLapse or
              EnabledEquations ==
                  Equations::HamiltonianLapseAndShift,
          tmpl::list<gr::Tags::Conformal<gr::Tags::StressTrace<DataVector>,
                                         ConformalMatterScale>,
                     ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataVector>>>,
          tmpl::list<>>,
      tmpl::conditional_t<
          EnabledEquations ==
              Equations::HamiltonianAndLapse,
          tmpl::list<
              Tags::LongitudinalShiftMinusDtConformalMetricSquare<DataVector>,
              Tags::ShiftDotDerivExtrinsicCurvatureTrace<DataVector>>,
          tmpl::list<>>,
      // Additional quantities for momentum constraint
      tmpl::conditional_t<
          EnabledEquations ==
              Equations::HamiltonianLapseAndShift,
          tmpl::list<
              gr::Tags::Conformal<gr::Tags::MomentumDensity<DataVector, 3>,
                                  ConformalMatterScale>,
              ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataVector>,
                            tmpl::size_t<3>, Frame::Inertial>,
              Tags::ShiftBackground<DataVector, 3, Frame::Inertial>,
              Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
                  DataVector, 3, Frame::Inertial>,
              // Note that this is the plain divergence, i.e. with no
              // Christoffel symbol terms added
              ::Tags::div<
                  Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
                      DataVector, 3, Frame::Inertial>>,
              tmpl::conditional_t<
                  ConformalGeometry == Geometry::Curved,
                  tmpl::list<
                      Tags::ConformalMetric<DataVector, 3, Frame::Inertial>,
                      Tags::ConformalChristoffelFirstKind<DataVector, 3,
                                                          Frame::Inertial>,
                      Tags::ConformalChristoffelSecondKind<DataVector, 3,
                                                           Frame::Inertial>>,
                  tmpl::list<>>>,
          tmpl::list<>>>>;
  using inv_metric_tag = tmpl::conditional_t<
      ConformalGeometry == Geometry::FlatCartesian, void,
      Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>>;

  using fluxes_computer = Fluxes<EnabledEquations, ConformalGeometry>;
  using sources_computer =
      Sources<EnabledEquations, ConformalGeometry, ConformalMatterScale>;
  using sources_computer_linearized =
      LinearizedSources<EnabledEquations, ConformalGeometry,
                        ConformalMatterScale>;

  using boundary_conditions_base =
      elliptic::BoundaryConditions::BoundaryCondition<3>;
};

}  // namespace Xcts
