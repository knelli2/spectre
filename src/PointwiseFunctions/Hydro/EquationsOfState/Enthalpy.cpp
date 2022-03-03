// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/Hydro/EquationsOfState/Enthalpy.hpp"

#include <cmath>
#include <functional>
#include <numeric>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/RootFinding/NewtonRaphson.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace EquationsOfState {
// Given an expansion h(x) = sum_i f_i(x) , compute int_a^x sum_i f_i(x) e^x +
// F(a)
Enthalpy::Coefficients Enthalpy::Coefficients::compute_exponential_integral(
    std::pair<double, double> initial_condition) {
  // Not yet implemented
  Enthalpy::Coefficients exponential_integral_coefficients = *this;
  return exponential_integral_coefficients;
}
Enthalpy::Coefficients Enthalpy::Coefficients::compute_derivative() {
  // Not yet implemeted
  Enthalpy::Coefficients derivative_coefficients = *this;
  return derivative_coefficients;
}
void Enthalpy::Coefficients::pup(PUP::er& p) {
  p | polynomial_coefficients;
  p | sin_coefficients;
  p | cos_coefficients;
  p | reference_density;
  p | trig_scale;
}

Enthalpy::Enthalpy(const double reference_density, const double max_density,
                   const double min_density, const double min_energy_density,
                   const double trig_scale,
                   const std::vector<double> polynomial_coefficients,
                   const std::vector<double> sin_coefficients,
                   const std::vector<double> cos_coefficients,
                   double lower_spectral_reference_density,
                   double lower_spectral_reference_pressure,
                   std::vector<double> lower_spectral_coefficients,
                   double lower_spectral_upper_density)
    : reference_density_(reference_density),
      minimum_density_(min_density),
      maximum_density_(max_density),
      lower_spectral_(
          lower_spectral_reference_density, lower_spectral_reference_pressure,
          lower_spectral_coefficients, lower_spectral_upper_density),
      coefficients_(polynomial_coefficients, sin_coefficients, cos_coefficients,
                    trig_scale, reference_density)

{
  minimum_enthalpy_ = specific_enthalpy_from_density(minimum_density_);
  exponential_integral_coefficients_ =
      coefficients_.compute_exponential_integral(
          {min_density, min_energy_density});
  derivative_coefficients_ = coefficients_.compute_derivative();
}

EQUATION_OF_STATE_MEMBER_DEFINITIONS(, Enthalpy, double, 1)
EQUATION_OF_STATE_MEMBER_DEFINITIONS(, Enthalpy, DataVector, 1)

Enthalpy::Enthalpy(CkMigrateMessage* /*unused*/) {}

void Enthalpy::pup(PUP::er& p) {
  EquationOfState<true, 1>::pup(p);
  p | reference_density_;
  p | lower_spectral_;
  p | maximum_density_;
  p | minimum_density_;
  p | coefficients_;
  p | exponential_integral_coefficients_;
  p | derivative_coefficients_;
}

double Enthalpy::x_from_density(const double density) const {
  return log(density / reference_density_);
};
double Enthalpy::density_from_x(const double x) const {
  return reference_density_ * exp(x);
};

double Enthalpy::energy_density_from_log_density(const double x) const {
  reference_density_* evaluate_coefficients(exponential_integral_coefficients_,
                                            x);
}

// Evaluate the function represented by the coefficinets at x  = log(rho)
double Enthalpy::evaluate_coefficients(
    const Enthalpy::Coefficients& coefficients, const double x) {
  // Not as good as Horner's rule, so don't use for polynomials
  auto evaluate_with_function_basis =
      +[](const std::function<double(double, int)>& basis_function,
          double evaluation_point,
          const std::vector<double>& basis_coefficients, double initial) {
        double sum = initial;
        for (size_t index = 0; index < basis_coefficients.size(); index++) {
          sum += basis_coefficients[index] *
                 basis_function(evaluation_point, index);
        }
        return sum;
      };
  const double polynomial_contribution =
      evaluate_polynomial(coefficients.polynomial_coefficients, x);
  const double sin_contribution = evaluate_with_function_basis(
      +[](double x, int n) { return sin((n + 1) * x); }, x,
      coefficients.sin_coefficients, 0.0);
  const double cos_contribution = evaluate_with_function_basis(
      +[](double x, int n) { return cos((n + 1) * x); }, x,
      coefficients.sin_coefficients, 0.0);
  double value =
      polynomial_contribution + (sin_contribution + cos_contribution);
  if (coefficients.exponential_prefactor) {
    value *= exp(x);
    value += coefficients.exponential_constant;
  }
  return value;
};

template <typename DataType>
Scalar<DataType> Spectral::pressure_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  if constexpr (std::is_same_v<DataType, double>) {
    return Scalar<double>{pressure_from_density(get(rest_mass_density))};
  } else if constexpr (std::is_same_v<DataType, DataVector>) {
    auto result = make_with_value<Scalar<DataVector>>(rest_mass_density, 0.0);
    for (size_t i = 0; i < get(result).size(); ++i) {
      get(result)[i] = pressure_from_density(get(rest_mass_density)[i]);
    }
    return result;
  }
}

template <class DataType>
Scalar<DataType> Spectral::rest_mass_density_from_enthalpy_impl(
    const Scalar<DataType>& specific_enthalpy) const {
  if constexpr (std::is_same_v<DataType, double>) {
    return Scalar<double>{
        rest_mass_density_from_enthalpy(get(specific_enthalpy))};
  } else if constexpr (std::is_same_v<DataType, DataVector>) {
    auto result = make_with_value<Scalar<DataVector>>(specific_enthalpy, 0.0);
    for (size_t i = 0; i < get(result).size(); ++i) {
      get(result)[i] =
          rest_mass_density_from_enthalpy(get(specific_enthalpy)[i]);
    }
    return result;
  }
}

template <class DataType>
Scalar<DataType> Spectral::specific_enthalpy_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  if constexpr (std::is_same_v<DataType, double>) {
    return Scalar<double>{
        specific_enthalpy_from_density(get(rest_mass_density))};
  } else if constexpr (std::is_same_v<DataType, DataVector>) {
    auto result = make_with_value<Scalar<DataVector>>(rest_mass_density, 0.0);
    for (size_t i = 0; i < get(result).size(); ++i) {
      get(result)[i] =
          specific_enthalpy_from_density(get(rest_mass_density)[i]);
    }
    return result;
  }
}

template <class DataType>
Scalar<DataType> Spectral::specific_internal_energy_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  if constexpr (std::is_same_v<DataType, double>) {
    return Scalar<double>{
        specific_internal_energy_from_density(get(rest_mass_density))};
  } else if constexpr (std::is_same_v<DataType, DataVector>) {
    auto result = make_with_value<Scalar<DataVector>>(rest_mass_density, 0.0);
    for (size_t i = 0; i < get(result).size(); ++i) {
      get(result)[i] =
          specific_internal_energy_from_density(get(rest_mass_density)[i]);
    }
    return result;
  }
}

template <class DataType>
Scalar<DataType> Spectral::chi_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  if constexpr (std::is_same_v<DataType, double>) {
    return Scalar<double>{chi_from_density(get(rest_mass_density))};
  } else if constexpr (std::is_same_v<DataType, DataVector>) {
    auto result = make_with_value<Scalar<DataVector>>(rest_mass_density, 0.0);
    for (size_t i = 0; i < get(result).size(); ++i) {
      get(result)[i] = chi_from_density(get(rest_mass_density)[i]);
    }
    return result;
  }
}

template <class DataType>
Scalar<DataType> Spectral::kappa_times_p_over_rho_squared_from_density_impl(
    const Scalar<DataType>& rest_mass_density) const {
  return make_with_value<Scalar<DataType>>(get(rest_mass_density), 0.0);
}

double Enthalpy::chi_from_density(const double rest_mass_density) const {
  const double x = log(rest_mass_density / reference_density_);
  return evaluate_coefficients(derivative_coefficients_, x);
}

double Enthalpy::specific_internal_energy_from_density(
    const double rest_mass_density) const {
  const double x = log(rest_mass_density / reference_density_);
  if (Enthalpy::in_spectral_domain(rest_mass_density)) {
    return lower_spectral_.specific_internal_energy_from_density(
        rest_mass_density);
  } else {
    return energy_density_from_log_density(x_from_density(rest_mass_density)) -
           rest_mass_density;
  }
}

double Enthalpy::specific_enthalpy_from_density(
    const double rest_mass_density) const {
  return evaluate_coefficients(coefficients_,
                               x_from_density(rest_mass_density));
}

// P(x) = rho(x)h(x) - e(x) with h the specific enthalpy, and e
// the energy density.
double Enthalpy::pressure_from_log_density(const double x) const {
  if (Enthalpy::in_spectral_domain(density_from_x(x))) {
    return lower_spectral_.pressure_from_log_density(x);
  } else {
    return density_from_x(x) * evaluate_coefficients(coefficients_, x) -
           energy_density_from_log_density(x);
  }
}

double Enthalpy::pressure_from_density(const double rest_mass_density) const {
  if (Enthalpy::in_spectral_domain(rest_mass_density)) {
    return lower_spectral_.pressure_from_density(rest_mass_density);
  } else {
    const double x = log(rest_mass_density / reference_density_);
    return pressure_from_log_density(x);
  }
}

// Solve for h(rho)=h0, which requires rootfinding for this EoS
double Enthalpy::rest_mass_density_from_enthalpy(
    const double specific_enthalpy) const {
  if (specific_enthalpy <= minimum_enthalpy_) {
    return lower_spectral_.rest_mass_density_from_enthalpy(specific_enthalpy);
  } else {
    // Root-finding appropriate between reference density and maximum density
    // We can use x=0 and x=x_max as bounds
    const auto f_df_lambda = [this, &specific_enthalpy](const double density) {
      const auto x = x_from_density(density);
      const double f =
          evaluate_coefficients(coefficients_, x) - specific_enthalpy;

      const double df = 1 / evaluate_coefficients(derivative_coefficients_, x);
      ;
      return std::make_pair(f, df);
    };
    const size_t digits = 14;
    const double intial_guess = 0.5 * (minimum_density_ + maximum_density_);
    const auto root_from_lambda = RootFinder::newton_raphson(
        f_df_lambda, intial_guess, reference_density_, maximum_density_,
        digits);
    return root_from_lambda;
  }
}

PUP::able::PUP_ID EquationsOfState::Enthalpy::my_PUP_ID = 0;

}  // namespace EquationsOfState
