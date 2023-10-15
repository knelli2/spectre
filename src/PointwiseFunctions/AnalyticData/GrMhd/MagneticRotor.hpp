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
 * \brief Analytic initial data for a magnetic rotor.
 *
 * This is a test first described in \cite Balsara1999 for classical MHD and
 * later generalised to relativistic MHD in \cite DelZanna2002rv
 *
 * This effectively 2D test initially consists of an infinitely long cylinder of
 * radius `RotorRadius` rotating about the z-axis with the given
 * `AngularVelocity`. The rest mass density of the fluid inside the rotor,
 * `RotorDensity`, is higher than the `BackgroundDensity` outside of the rotor.
 * The fluid is at a constant `Pressure`.  The rotor is embedded in a constant
 * `MagneticField` (usually taken to be along the x-axis).  The fluid is an
 * ideal fluid with the given `AdiabaticIndex`.  Evolving the initial data,
 * magnetic braking will slow down the rotor, while dragging the magnetic field
 * lines.
 *
 * The standard test setup is done on a unit cube \f$[-0.5,0.5]^3\f$ with the
 * following values given for the options:
 * -  RotorRadius: 0.1
 * -  RotorDensity: 10.0
 * -  BackgroundDensity: 1.0
 * -  Pressure: 1.0
 * -  AngularVelocity: 9.95
 * -  MagneticField: [1.0, 0.0, 0.0]
 * -  AdiabaticIndex: 1.66666666666666667
 *
 * Note that \cite Zanotti2016efficient uses different parameters,
 * -  RotorRadius: 0.1
 * -  RotorDensity: 10.0
 * -  BackgroundDensity: 1.0
 * -  Pressure: 1.0
 * -  AngularVelocity: 9.3
 * -  MagneticField: [1.0, 0.0, 0.0]
 * -  AdiabaticIndex: 1.333333333333333
 *
 * The magnetic field in the disk should rotate by about 90 degrees by t = 0.4.
 */
class MagneticRotor : public evolution::initial_data::InitialData,
                      public MarkAsAnalyticData,
                      public AnalyticDataBase,
                      public hydro::TemperatureInitialization<MagneticRotor> {
 public:
  using equation_of_state_type = EquationsOfState::IdealFluid<true>;

  /// Radius of the rotor.
  struct RotorRadius {
    using type = double;
    inline const static std::string help {
        "The initial radius of the rotor."};
    static type lower_bound() { return 0.0; }
  };
  /// Density inside the rotor.
  struct RotorDensity {
    using type = double;
    inline const static std::string help {"Density inside RotorRadius."};
    static type lower_bound() { return 0.0; }
  };
  /// Density outside the rotor.
  struct BackgroundDensity {
    using type = double;
    inline const static std::string help {"Density outside RotorRadius."};
    static type lower_bound() { return 0.0; }
  };
  /// Uniform pressure inside and outside the rotor.
  struct Pressure {
    using type = double;
    inline const static std::string help {"Pressure."};
    static type lower_bound() { return 0.0; }
  };
  /// Angular velocity inside the rotor.
  struct AngularVelocity {
    using type = double;
    inline const static std::string help {
        "Angular velocity of matter inside RotorRadius"};
  };
  /// The x,y,z components of the uniform magnetic field threading the matter.
  struct MagneticField {
    using type = std::array<double, 3>;
    inline const static std::string help {
        "The x,y,z components of the uniform magnetic field."};
  };
  /// The adiabatic index of the ideal fluid.
  struct AdiabaticIndex {
    using type = double;
    inline const static std::string help {
        "The adiabatic index of the ideal fluid."};
    static type lower_bound() { return 1.0; }
  };

  using options =
      tmpl::list<RotorRadius, RotorDensity, BackgroundDensity, Pressure,
                 AngularVelocity, MagneticField, AdiabaticIndex>;

  inline const static std::string help {
      "Magnetic rotor analytic initial data."};

  MagneticRotor() = default;
  MagneticRotor(const MagneticRotor& /*rhs*/) = default;
  MagneticRotor& operator=(const MagneticRotor& /*rhs*/) = default;
  MagneticRotor(MagneticRotor&& /*rhs*/) = default;
  MagneticRotor& operator=(MagneticRotor&& /*rhs*/) = default;
  ~MagneticRotor() override = default;

  MagneticRotor(double rotor_radius, double rotor_density,
                double background_density, double pressure,
                double angular_velocity,
                const std::array<double, 3>& magnetic_field,
                double adiabatic_index, const Options::Context& context = {});

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit MagneticRotor(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(MagneticRotor);
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
  double rotor_radius_ = std::numeric_limits<double>::signaling_NaN();
  double rotor_density_ = std::numeric_limits<double>::signaling_NaN();
  double background_density_ = std::numeric_limits<double>::signaling_NaN();
  double pressure_ = std::numeric_limits<double>::signaling_NaN();
  double angular_velocity_ = std::numeric_limits<double>::signaling_NaN();
  std::array<double, 3> magnetic_field_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  double adiabatic_index_ = std::numeric_limits<double>::signaling_NaN();
  EquationsOfState::IdealFluid<true> equation_of_state_{};
  gr::Solutions::Minkowski<3> background_spacetime_{};

  friend bool operator==(const MagneticRotor& lhs, const MagneticRotor& rhs);

  friend bool operator!=(const MagneticRotor& lhs, const MagneticRotor& rhs);
};

}  // namespace grmhd::AnalyticData
