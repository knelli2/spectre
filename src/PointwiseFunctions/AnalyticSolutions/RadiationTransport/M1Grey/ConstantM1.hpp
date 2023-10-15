// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <pup.h>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"  // IWYU pragma: keep
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"  // IWYU pragma: keep
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma: no_forward_declare Tensor

/// \cond
class DataVector;
/// \endcond

namespace RadiationTransport {
namespace M1Grey {
namespace Solutions {

/*!
 * \brief Constant solution to M1 equations in Minkowski spacetime.
 *
 * An analytic solution to the 3-D M1 system. The user specifies the mean
 * flow velocity of the fluid and the radiation energy density in the
 * fluid frame J.
 * The radiation is taken to be in equilibrium with the fluid
 * (i.e. comoving, with an isotropic pressure P=J/3)
 *
 */
class ConstantM1 : public MarkAsAnalyticSolution {
 public:
  /// The mean flow velocity.
  struct MeanVelocity {
    using type = std::array<double, 3>;
    inline const static std::string help {"The mean flow velocity."};
  };

  /// The radiation comoving energy density
  struct ComovingEnergyDensity {
    using type = double;
    inline const static std::string help {
        "Comoving energy density of radiation."};
  };

  using options = tmpl::list<MeanVelocity, ComovingEnergyDensity>;
  inline const static std::string help {
      "Constant radiation field in equilibrium with fluid, in Minkowski "
      "spacetime."};

  ConstantM1() = default;
  ConstantM1(const ConstantM1& /*rhs*/) = delete;
  ConstantM1& operator=(const ConstantM1& /*rhs*/) = delete;
  ConstantM1(ConstantM1&& /*rhs*/) = default;
  ConstantM1& operator=(ConstantM1&& /*rhs*/) = default;
  ~ConstantM1() = default;

  ConstantM1(const std::array<double, 3>& mean_velocity,
             double comoving_energy_density);

  explicit ConstantM1(CkMigrateMessage* /*unused*/) {}

  /// @{
  /// Retrieve fluid and neutrino variables at `(x, t)`
  template <typename NeutrinoSpecies>
  auto variables(const tnsr::I<DataVector, 3>& x, double t,
                 tmpl::list<RadiationTransport::M1Grey::Tags::TildeE<
                     Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeE<
          Frame::Inertial, NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(const tnsr::I<DataVector, 3>& x, double t,
                 tmpl::list<RadiationTransport::M1Grey::Tags::TildeS<
                     Frame::Inertial, NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::TildeS<
          Frame::Inertial, NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(const tnsr::I<DataVector, 3>& x, double t,
                 tmpl::list<RadiationTransport::M1Grey::Tags::GreyEmissivity<
                     NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<
          RadiationTransport::M1Grey::Tags::GreyEmissivity<NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(
      const tnsr::I<DataVector, 3>& x, double t,
      tmpl::list<RadiationTransport::M1Grey::Tags::GreyAbsorptionOpacity<
          NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                                 GreyAbsorptionOpacity<NeutrinoSpecies>>;

  template <typename NeutrinoSpecies>
  auto variables(
      const tnsr::I<DataVector, 3>& x, double t,
      tmpl::list<RadiationTransport::M1Grey::Tags::GreyScatteringOpacity<
          NeutrinoSpecies>> /*meta*/) const
      -> tuples::TaggedTuple<RadiationTransport::M1Grey::Tags::
                                 GreyScatteringOpacity<NeutrinoSpecies>>;

  auto variables(const tnsr::I<DataVector, 3>& x, double t,
                 tmpl::list<hydro::Tags::LorentzFactor<DataVector>> /*meta*/)
      const -> tuples::TaggedTuple<hydro::Tags::LorentzFactor<DataVector>>;

  auto variables(
      const tnsr::I<DataVector, 3>& x, double t,
      tmpl::list<hydro::Tags::SpatialVelocity<DataVector, 3>> /*meta*/) const
      -> tuples::TaggedTuple<hydro::Tags::SpatialVelocity<DataVector, 3>>;
  /// @}

  /// Retrieve a collection of fluid and neutrino variables at `(x, t)`
  template <typename... Tags>
  tuples::TaggedTuple<Tags...> variables(const tnsr::I<DataVector, 3>& x,
                                         double t,
                                         tmpl::list<Tags...> /*meta*/) const {
    static_assert(sizeof...(Tags) > 1,
                  "The generic template will recurse infinitely if only one "
                  "tag is being retrieved.");
    return {get<Tags>(variables(x, t, tmpl::list<Tags>{}))...};
  }

  /// Retrieve the metric variables
  template <typename Tag>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataVector, 3>& x, double t,
                                     tmpl::list<Tag> /*meta*/) const {
    return background_spacetime_.variables(x, t, tmpl::list<Tag>{});
  }

  // clang-tidy: no runtime references
  void pup(PUP::er& /*p*/);  //  NOLINT

 private:
  friend bool operator==(const ConstantM1& lhs, const ConstantM1& rhs);

  std::array<double, 3> mean_velocity_{
      {std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN(),
       std::numeric_limits<double>::signaling_NaN()}};
  double comoving_energy_density_ =
      std::numeric_limits<double>::signaling_NaN();
  gr::Solutions::Minkowski<3> background_spacetime_{};
};

bool operator!=(const ConstantM1& lhs, const ConstantM1& rhs);
}  // namespace Solutions
}  // namespace M1Grey
}  // namespace RadiationTransport
