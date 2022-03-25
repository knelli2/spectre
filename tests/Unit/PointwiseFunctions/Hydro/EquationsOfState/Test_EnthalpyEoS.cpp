// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <limits>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestCreation.hpp"
#include "Helpers/PointwiseFunctions/Hydro/EquationsOfState/TestHelpers.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/Enthalpy.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"

namespace {

void check_exact() {
  // Test creation
  // TODO : Find an analytic EoS that isn't weird for test case
  namespace EoS = EquationsOfState;
  Parallel::printf("Failed after entry");
  const std::vector<double> poly_coefs = {1.0, 0.2, 0.0, 0.0, 0.0001};
  const std::vector<double> sin_coefs = {0.01, 0.003, -0.0001, .0001};
  const std::vector<double> cos_coefs = {0.01, 0.003, 0.0001, 0.00001};
  const double min_density = 4.0;
  const double max_density = 100.0;
  const double trig_scaling = 1.5;
  const double reference_density = 2.0;
  const double min_energy_density = 4.448303974107519;
  const double lower_spectral_reference_density = .8;
  const double lower_spectral_reference_pressure = 0.3727074036491289;
  const double lower_spectral_upper_density = 4;
  const std::vector<double> lower_spectral_gamma_coefficients{.30255164};
  EoS::Spectral lower_spectral{
      lower_spectral_reference_density, lower_spectral_reference_pressure,
      lower_spectral_gamma_coefficients, lower_spectral_upper_density};
  Parallel::printf("Made the Spectral but not the Enthalpy");
  // Parallel::register_derived_classes_with_charm<
  //   EoS::EquationOfState<true, 1>>();
  TestHelpers::test_creation<std::unique_ptr<EoS::EquationOfState<true, 1>>>(
      {"Enthalpy:\n"
       "  ReferenceDensity: 2.0\n"
       "  MaximumDensity: 4.0\n"
       "  MinimumDensity: 100.0\n"
       "  MinimumEnergyDensity: 4.448303974107519\n"
       "  TrigScaling: 1.5\n"
       "  PolynomialCoefficients: [1.0,0.2,0.0,0.0,0.0001]\n"
       "  SinCoefficients: [0.01,0.003,-0.0001,0.0001]\n"
       "  CosCoefficients: [0.01,0.003,0.0001,0.00001]\n"
       "  Spectral:\n"
       "    ReferenceDensity: 0.8\n"
       "    ReferencePressure: 0.3727074036491289\n"
       "    Coefficients: [0.30255164]\n"
       "    UpperDensity: 4.0\n"});
  Parallel::printf("Initiated test thing");
  EoS::Enthalpy<EoS::Spectral> eos(reference_density, max_density, min_density,
                                   min_energy_density, trig_scaling, poly_coefs,
                                   sin_coefs, cos_coefs, lower_spectral);
  // Test DataVector functions
  {
    Parallel::printf("Failed after lots of stuff but before any tests");
    const Scalar<DataVector> rho{DataVector{1.5 * exp(1.0), 1.5 * exp(2.0),
                                            1.5 * exp(3.0), 1.5 * exp(4.0)}};
    const Scalar<DataVector> p = eos.pressure_from_density(rho);
    INFO(rho);
    INFO(p);
    // TODO, put the actual right values in here
    const Scalar<DataVector> p_expected{
        DataVector{0.1790434090653277, 1.46117617057407734, 5.2223018314642724,
                   17.05980881997609089}};
    CHECK_ITERABLE_APPROX(p, p_expected);
    const auto eps_c = eos.specific_internal_energy_from_density(rho);
    const Scalar<DataVector> eps_expected{
        DataVector{0.11289220759026733, 0.20684752054286823,
                   0.36264928982309969, 0.55159228289061324}};
    CHECK_ITERABLE_APPROX(eps_c, eps_expected);
    const auto h_c = eos.specific_enthalpy_from_density(rho);
    const auto h_expected =
        get(eps_expected) + 1.0 + get(p_expected) / get(rho);
    CHECK_ITERABLE_APPROX(get(h_c), h_expected);
    const auto chi_c = eos.chi_from_density(rho);
    const Scalar<DataVector> chi_expected{
        DataVector{0.18207356789438428, 0.1923165022214991, 0.19908149093638267,
                   0.25186554704031844}};
    CHECK_ITERABLE_APPROX(get(chi_expected), get(chi_c));
    const Scalar<DataVector> p_c_kappa_c_over_rho_sq_expected{
        DataVector{0.0, 0.0, 0.0, 0.0}};
    const auto p_c_kappa_c_over_rho_sq =
        eos.kappa_times_p_over_rho_squared_from_density(rho);
    CHECK_ITERABLE_APPROX(p_c_kappa_c_over_rho_sq_expected,
                          p_c_kappa_c_over_rho_sq);
    const auto rho_from_enthalpy = eos.rest_mass_density_from_enthalpy(h_c);
    CHECK_ITERABLE_APPROX(rho, rho_from_enthalpy);
  }
  // Test double functions
  {
    const Scalar<double> rho{1.5 * exp(1.0)};
    const auto p = eos.pressure_from_density(rho);
    const double p_expected = 0.1790434090653277;
    CHECK(get(p) == approx(p_expected));
    const auto eps = eos.specific_internal_energy_from_density(rho);
    const double eps_expected = 0.11289220759026733;
    CHECK_ITERABLE_APPROX(get(eps), eps_expected);
    const auto h = eos.specific_enthalpy_from_density(rho);
    const double h_expected = eps_expected + 1.0 + p_expected / get(rho);
    CHECK_ITERABLE_APPROX(get(h), h_expected);
    const auto chi = eos.chi_from_density(rho);
    const double chi_expected = 0.18207356789438428;
    CHECK_ITERABLE_APPROX(chi_expected, get(chi));
    const auto p_c_kappa_c_over_rho_sq =
        eos.kappa_times_p_over_rho_squared_from_density(rho);
    const double p_c_kappa_c_over_rho_sq_expected = 0.0;
    CHECK_ITERABLE_APPROX(p_c_kappa_c_over_rho_sq_expected,
                          get(p_c_kappa_c_over_rho_sq));
    const auto rho_from_enthalpy = eos.rest_mass_density_from_enthalpy(h);
    CHECK_ITERABLE_APPROX(get(rho), get(rho_from_enthalpy));
  }
  // Test bounds
  CHECK(0.0 == eos.rest_mass_density_lower_bound());
  CHECK(0.0 == eos.specific_internal_energy_lower_bound(1.0));
  CHECK(1.0 == eos.specific_enthalpy_lower_bound());
  const double max_double = std::numeric_limits<double>::max();
  CHECK(max_double == eos.rest_mass_density_upper_bound());
  CHECK(max_double == eos.specific_internal_energy_upper_bound(1.0));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.EquationsOfState.Enthalpy",
                  "[Unit][EquationsOfState]") {
  check_exact();
}
