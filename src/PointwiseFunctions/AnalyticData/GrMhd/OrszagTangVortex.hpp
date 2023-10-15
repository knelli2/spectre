// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/IdealFluid.hpp"  // IWYU pragma: keep
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma: no_include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"  // for IdealFluid

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace grmhd::AnalyticData {

/*!
 * \brief Analytic initial data for the relativistic Orszag-Tang vortex.
 *
 * The relativistic version of the Orszag-Tang vortex is a
 * 2-dimensional test case for relativistic MHD systems (see, e.g.,
 * \cite Beckwith2011iy).  It describes the flow of an ideal fluid with
 * adiabatic index \f$5/3\f$.  The initial conditions (and hence the
 * states at later times) are periodic in both \f$x\f$ and \f$y\f$
 * with period 1.  The initial conditions are:
 * \f{align*}
 * \rho &= \frac{25}{36 \pi} \\
 * p &= \frac{5}{12 \pi} \\
 * v_x &= -\frac{1}{2} \sin(2 \pi y) \\
 * v_y &= \frac{1}{2} \sin(2 \pi x) \\
 * B_x &= -\frac{1}{\sqrt{4 \pi}} \sin(2 \pi y) \\
 * B_y &= \frac{1}{\sqrt{4 \pi}} \sin(4 \pi x)
 * \f}
 * with \f$\rho\f$ the rest mass density, \f$p\f$ the pressure,
 * \f$v_i\f$ the spatial velocity, and \f$B_i\f$ the magnetic field.
 *
 * \parblock
 * \note We do not currently support 2-dimensional RMHD, so this class
 * provides 3-dimensional data with no \f$z\f$-dependence.
 * \endparblock
 *
 * \parblock
 * \note There are multiple errors in the description of this test
 * problem in the original SpECTRE paper \cite Kidder2016hev and there
 * is a sign error in the velocity in \cite Beckwith2011iy.  Despite these
 * errors, the actual tests performed for those papers matched the standard
 * problem as presented here.
 * \endparblock
 */
class OrszagTangVortex
    : public evolution::initial_data::InitialData,
      public AnalyticDataBase,
      public hydro::TemperatureInitialization<OrszagTangVortex>,
      public MarkAsAnalyticData {
 public:
  using equation_of_state_type = EquationsOfState::IdealFluid<true>;

  using options = tmpl::list<>;

  inline const static std::string help {
      "The relativistic Orszag-Tang vortex"};

  OrszagTangVortex();
  OrszagTangVortex(const OrszagTangVortex& /*rhs*/) = default;
  OrszagTangVortex& operator=(const OrszagTangVortex& /*rhs*/) = default;
  OrszagTangVortex(OrszagTangVortex&& /*rhs*/) = default;
  OrszagTangVortex& operator=(OrszagTangVortex&& /*rhs*/) = default;
  ~OrszagTangVortex() override = default;

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit OrszagTangVortex(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(OrszagTangVortex);
  /// \endcond

  /// @{
  /// Retrieve hydro variable at `x`
  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::RestMassDensity<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::RestMassDensity<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::ElectronFraction<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::ElectronFraction<DataType>>;

  template <typename DataType>
  auto variables(
      const tnsr::I<DataType, 3>& x,
      tmpl::list<hydro::Tags::SpecificInternalEnergy<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::SpecificInternalEnergy<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::Pressure<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::Pressure<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::SpatialVelocity<DataType, 3>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataType, 3>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::MagneticField<DataType, 3>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::MagneticField<DataType, 3>>;

  template <typename DataType>
  auto variables(
      const tnsr::I<DataType, 3>& x,
      tmpl::list<hydro::Tags::DivergenceCleaningField<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::DivergenceCleaningField<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::LorentzFactor<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::SpecificEnthalpy<DataType>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::SpecificEnthalpy<DataType>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x,
                 tmpl::list<hydro::Tags::Temperature<DataType>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::Temperature<DataType>> {
    return TemperatureInitialization::variables(
        x, tmpl::list<hydro::Tags::Temperature<DataType>>{});
  }
  /// @}

  /// Retrieve a collection of hydrodynamic variables at position x
  template <typename DataType, typename Tag1, typename Tag2, typename... Tags>
  tuples::TaggedTuple<Tag1, Tag2, Tags...> variables(
      const tnsr::I<DataType, 3>& x,
      tmpl::list<Tag1, Tag2, Tags...> /*meta*/) const {
    return {tuples::get<Tag1>(variables(x, tmpl::list<Tag1>{})),
            tuples::get<Tag2>(variables(x, tmpl::list<Tag2>{})),
            tuples::get<Tags>(variables(x, tmpl::list<Tags>{}))...};
  }

  /// Retrieve the metric variables
  template <typename DataType, typename Tag,
            Requires<not tmpl::list_contains_v<hydro::grmhd_tags<DataType>,
                                               Tag>> = nullptr>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataType, 3>& x,
                                     tmpl::list<Tag> /*meta*/) const {
    constexpr double dummy_time = 0.;
    return {std::move(get<Tag>(gr::Solutions::Minkowski<3>{}.variables(
        x, dummy_time, tmpl::list<Tag>{})))};
  }

  const equation_of_state_type& equation_of_state() const {
    return equation_of_state_;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override;

 private:
  EquationsOfState::IdealFluid<true> equation_of_state_{5. / 3.};
};

bool operator==(const OrszagTangVortex& lhs, const OrszagTangVortex& rhs);

bool operator!=(const OrszagTangVortex& lhs, const OrszagTangVortex& rhs);

}  // namespace grmhd::AnalyticData
