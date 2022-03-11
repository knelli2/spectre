// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/list/adt.hpp>
#include <boost/preprocessor/repetition/for.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/tuple/to_list.hpp>
#include <limits>
#include <pup.h>
#include <vector>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"  // IWYU pragma: keep
#include "PointwiseFunctions/Hydro/EquationsOfState/Spectral.hpp"
#include "Utilities/Math.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;

/// \endcond

// IWYU pragma: no_forward_declare Tensor

namespace EquationsOfState {

// Convenience from not having to change all of the code everytime the details
// of the implementation change

/*!
 * \ingroup EquationsOfStateGroup
 * \brief An equation of state given by parametrized enthalpy
 *
 * This equation of state is determined as a function of \f$x =
 * \ln(\rho/\rho_0)\f$ where \f$\rho\f$ is the rest mass density and
 * \f$\rho_0\f$ is the provided reference density.
 * The pseudo-enthalpy \f$h \equiv (p + rho  + epsilon)/rho\f$
 * is expaded as
 *
 * \f{equation}
 * \h(x) = \sum_i a_i x^i + \sum_j b_j \sin(jkx) + c_j \cos(jkx)
 * \f}
 * for the set of spectral coefficinets \f$\gamma_n\f$ when
 * \f$0 < x < x_u = \ln(\rho_u/\rho_0)\f$, where \f$\rho_u\f$ is the provided
 * upper density.
 *
 * For \f$ x < 0 \f$, \f$ \Gamma(x) = \gamma_0 \f$.
 *
 * For \f$ x > x_u \f$, \f$ \Gamma(x) = \Gamma(x_u) \f$
 *
 *
 */
class Enthalpy : public EquationOfState<true, 1> {
 private:
  struct Coefficients {
    std::vector<double> polynomial_coefficients;
    std::vector<double> sin_coefficients;
    std::vector<double> cos_coefficients;
    double trig_scale;
    double reference_density;
    bool has_exponential_prefactor;
    double exponential_external_constant;
    Coefficients() = default;
    Coefficients(const Coefficients& coefficients) = default;
    Coefficients(std::vector<double> in_polynomial_coefficients,
                 std::vector<double> in_sin_coefficients,
                 std::vector<double> in_cos_coefficients, double in_trig_scale,
                 double in_reference_density,
                 double in_exponential_constant =
                     std::numeric_limits<double>::quiet_NaN());

    Enthalpy::Coefficients compute_exponential_integral(
        std::pair<double, double> initial_condition);
    Enthalpy::Coefficients compute_derivative();
    void pup(PUP::er& p);
  };

 public:
  static constexpr size_t thermodynamic_dim = 1;
  static constexpr bool is_relativistic = true;

  struct ReferenceDensity {
    using type = double;
    static constexpr Options::String help = {"Reference density rho_0"};
    static double lower_bound() { return 0.0; }
  };

  struct MinimumDensity {
    using type = double;
    static constexpr Options::String help = {
        "Minimum valid density rho_min,"
        " for this parametrization"};
    static double lower_bound() { return 0.0; }
  };
  struct MaximumDensity {
    using type = double;
    static constexpr Options::String help = {"Upper density rho_u"};
    static double lower_bound() { return 0.0; }
  };

  struct PolynomialCoefficients {
    using type = std::vector<double>;
    static constexpr Options::String help = {"Polynomial coefficients a_i"};
  };

  struct TrigScaling {
    using type = double;
    static constexpr Options::String help = {
        "Fundamental wavenumber of trig "
        "functions, k"};
    static double lower_bound() { return 0.0; }
  };

  struct SinCoefficients {
    using type = std::vector<double>;
    static constexpr Options::String help = {"Sine coefficients b_j"};
  };
  struct CosCoefficients {
    using type = std::vector<double>;
    static constexpr Options::String help = {"Cosine coefficients c_j"};
  };
  struct LowDensitySpectral {
    using type = EquationsOfState::Spectral;
    static constexpr Options::String help = {
        "Low density EoS stitched at the MinimumDensity"};
  };

  static constexpr Options::String help = {
      "An EoS with a parametrized value h(log(rho/rho_0)) with h the specific "
      "enthalpy and rho the baryon rest mass density.  The enthalpy is "
      "expanded as a sum of polynomial terms and trigonometric corrections. "
      "let x = log(rho/rho_0) in"
      "h(x) = \\sum_i a_ix^i + \\sum_j b_jsin(k * j * x) + c_jcos(k * j * x) "
      "Note that rho(x)(1+epsilon(x)) = int_0^x e^x' h((x') dx' can be "
      "computed "
      "analytically, and therefore so can "
      "P(x) = rho(x) * (h(x) - (1 + epsilon(x))) "};

  using options = tmpl::list<ReferenceDensity, MaximumDensity, MinimumDensity,
                             PolynomialCoefficients, SinCoefficients,
                             CosCoefficients, LowDensitySpectral>;

  Enthalpy() = default;
  Enthalpy(const Enthalpy&) = default;
  Enthalpy& operator=(const Enthalpy&) = default;
  Enthalpy(Enthalpy&&) = default;
  Enthalpy& operator=(Enthalpy&&) = default;
  ~Enthalpy() override = default;

  Enthalpy(double reference_density, double max_density, double min_density,
           double min_energy_density, double trig_scale,
           std::vector<double> polynomial_coefficients,
           std::vector<double> sin_coefficients,
           std::vector<double> cos_coefficients,
           const EquationsOfState::Spectral& lower_spectral);

  EQUATION_OF_STATE_FORWARD_DECLARE_MEMBERS(Enthalpy, 1)

  WRAPPED_PUPable_decl_base_template(  // NOLINT
      SINGLE_ARG(EquationOfState<true, 1>), Enthalpy);

  /// The lower bound of the rest mass density that is valid for this EOS
  double rest_mass_density_lower_bound() const override { return 0.0; }

  /// The upper bound of the rest mass density that is valid for this EOS
  double rest_mass_density_upper_bound() const override {
    return std::numeric_limits<double>::max();
  }

  /// The lower bound of the specific internal energy that is valid for this EOS
  /// at the given rest mass density \f$\rho\f$
  double specific_internal_energy_lower_bound(
      const double /* rest_mass_density */) const override {
    return 0.0;
  }

  /// The upper bound of the specific internal energy that is valid for this EOS
  /// at the given rest mass density \f$\rho\f$
  double specific_internal_energy_upper_bound(
      const double /* rest_mass_density */) const override {
    return std::numeric_limits<double>::max();
  }

  /// The lower bound of the specific enthalpy that is valid for this EOS
  double specific_enthalpy_lower_bound() const override { return 1.0; }

 private:
  EQUATION_OF_STATE_FORWARD_DECLARE_MEMBER_IMPLS(1)

  bool in_spectral_domain(const double density) const {
    return density < minimum_density_;
  }

  double x_from_density(const double density) const;
  double density_from_x(const double x) const;
  double energy_density_from_log_density(const double x) const;

  static double evaluate_coefficients(
      const Enthalpy::Coefficients& coefficients, const double x);
  double chi_from_density(const double density) const;
  double specific_internal_energy_from_density(const double density) const;
  double specific_enthalpy_from_density(const double density) const;
  double pressure_from_density(const double density) const;
  double pressure_from_log_density(const double x) const;
  double rest_mass_density_from_enthalpy(const double specific_enthalpy) const;

  double reference_density_ = std::numeric_limits<double>::signaling_NaN();
  double minimum_density_ = std::numeric_limits<double>::signaling_NaN();
  double maximum_density_ = std::numeric_limits<double>::signaling_NaN();
  double minimum_enthalpy_ = std::numeric_limits<double>::signaling_NaN();

  EquationsOfState::Spectral lower_spectral_;
  Coefficients coefficients_;
  Coefficients exponential_integral_coefficients_;
  Coefficients derivative_coefficients_;
};

}  // namespace EquationsOfState
