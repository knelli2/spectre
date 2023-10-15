// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/IdealFluid.hpp"  // IWYU pragma: keep
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma:  no_include <pup.h>

/// \cond
namespace PUP {
class er;  // IWYU pragma: keep
}  // namespace PUP
/// \endcond

namespace grmhd::AnalyticData {

/*!
 * \brief Analytic initial data for an advecting magnetic field loop.
 *
 * This test, originally proposed in \cite Gardiner2005hy and presented in a
 * slightly modified form by \cite Mignone2010br, a region with annular
 * cross section with the specified `InnerRadius` and `OuterRadius` is given a
 * non-zero azimuthal magnetic field of constant magnitude `MagFieldStrength`
 * with zero magnetic field outside the loop.  Inside the `InnerRadius` the
 * magnetic field strength falls to zero quadratically. The loop is embedded in
 * an ideal fluid with the given `AdiabaticIndex`, `RestMassDensity` and
 * `Pressure` with a uniform `AdvectionVelocity`.  The magnetic field loop
 * should advect across the grid, maintaining its shape and strength, as long
 * as the magnetic pressure is negligible compared to the thermal pressure.
 *
 * This test diagnoses how well the evolution scheme preserves the no-monopole
 * condition, as well as the diffusivity of the scheme.
 *
 * The standard test setup is done on \f$x \in [-1,1]\f$, \f$y \in [-0.5,
 * 0.5]\f$, with periodic boundary conditions and with the following values
 * given for the options:
 * -  InnerRadius: 0.06
 * -  OuterRadius: 0.3
 * -  RestMassDensity: 1.0
 * -  Pressure: 1.0
 * -  AdvectionVelocity: [0.08164965809277261, 0.040824829046386304, 0.0]
 * -  MagFieldStrength: 0.001
 * -  AdiabaticIndex: 1.66666666666666667
 *
 */
class MagneticFieldLoop
    : public evolution::initial_data::InitialData,
      public MarkAsAnalyticData,
      public AnalyticDataBase,
      public hydro::TemperatureInitialization<MagneticFieldLoop> {
 public:
  using equation_of_state_type = EquationsOfState::IdealFluid<true>;

  /// The pressure throughout the fluid.
  struct Pressure {
    using type = double;
    inline const static std::string help {
        "The constant pressure throughout the fluid."};
    static type lower_bound() { return 0.0; }
  };

  /// The rest mass density throughout the fluid.
  struct RestMassDensity {
    using type = double;
    inline const static std::string help {
        "The constant density throughout the fluid."};
    static type lower_bound() { return 0.0; }
  };

  /// The adiabatic index for the ideal fluid.
  struct AdiabaticIndex {
    using type = double;
    inline const static std::string help {
        "The adiabatic index for the ideal fluid."};
    static type lower_bound() { return 1.0; }
  };

  /// The fluid velocity.
  struct AdvectionVelocity {
    using type = std::array<double, 3>;
    inline const static std::string help {"The advection velocity."};
    static type lower_bound() { return {{-1.0, -1.0, -1.0}}; }
    static type upper_bound() { return {{1.0, 1.0, 1.0}}; }
  };

  /// The strength of the magnetic field.
  struct MagFieldStrength {
    using type = double;
    inline const static std::string help {
        "The magnitude of the magnetic field."};
    static type lower_bound() { return 0.0; }
  };

  /// The inner radius of the magnetic loop.
  struct InnerRadius {
    using type = double;
    inline const static std::string help {
        "The inner radius of the magnetic loop."};
    static type lower_bound() { return 0.0; }
  };

  /// The outer radius of the magnetic loop.
  struct OuterRadius {
    using type = double;
    inline const static std::string help {
        "The outer radius of the magnetic loop."};
    static type lower_bound() { return 0.0; }
  };

  using options =
      tmpl::list<Pressure, RestMassDensity, AdiabaticIndex, AdvectionVelocity,
                 MagFieldStrength, InnerRadius, OuterRadius>;
  inline const static std::string help {
      "Periodic advection of a magnetic field loop in Minkowski."};

  MagneticFieldLoop() = default;
  MagneticFieldLoop(const MagneticFieldLoop& /*rhs*/) = default;
  MagneticFieldLoop& operator=(const MagneticFieldLoop& /*rhs*/) = default;
  MagneticFieldLoop(MagneticFieldLoop&& /*rhs*/) = default;
  MagneticFieldLoop& operator=(MagneticFieldLoop&& /*rhs*/) = default;
  ~MagneticFieldLoop() override = default;

  MagneticFieldLoop(double pressure, double rest_mass_density,
                    double adiabatic_index,
                    const std::array<double, 3>& advection_velocity,
                    double magnetic_field_magnitude, double inner_radius,
                    double outer_radius, const Options::Context& context = {});

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit MagneticFieldLoop(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(MagneticFieldLoop);
  /// \endcond

  /// @{
  /// Retrieve the GRMHD variables at a given position.
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
            Requires<tmpl::list_contains_v<
                gr::analytic_solution_tags<3, DataType>, Tag>> = nullptr>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataType, 3>& x,
                                     tmpl::list<Tag> /*meta*/) const {
    constexpr double dummy_time = 0.0;
    return background_spacetime_.variables(x, dummy_time, tmpl::list<Tag>{});
  }

  const EquationsOfState::IdealFluid<true>& equation_of_state() const {
    return equation_of_state_;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override;

 private:
  double pressure_ = std::numeric_limits<double>::signaling_NaN();
  double rest_mass_density_ = std::numeric_limits<double>::signaling_NaN();
  double adiabatic_index_ = std::numeric_limits<double>::signaling_NaN();
  std::array<double, 3> advection_velocity_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  double magnetic_field_magnitude_ =
      std::numeric_limits<double>::signaling_NaN();
  double inner_radius_ = std::numeric_limits<double>::signaling_NaN();
  double outer_radius_ = std::numeric_limits<double>::signaling_NaN();

  EquationsOfState::IdealFluid<true> equation_of_state_{};
  gr::Solutions::Minkowski<3> background_spacetime_{};

  friend bool operator==(const MagneticFieldLoop& lhs,
                         const MagneticFieldLoop& rhs);

  friend bool operator!=(const MagneticFieldLoop& lhs,
                         const MagneticFieldLoop& rhs);
};

}  // namespace grmhd::AnalyticData
