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

// namespace{
// // Is there an actual factorial function somewhere?
// size_t factorial(size_t n){
//   if (n == 0){return static_cast<size_t>(1);}
//   else{
//     return n * factorial(n - 1);
//   }
// }
// }

namespace EquationsOfState {
Enthalpy::Coefficients::Coefficients(
    std::vector<double> in_polynomial_coefficients,
    std::vector<double> in_sin_coefficients,
    std::vector<double> in_cos_coefficients, double in_trig_scale,
    double in_reference_density, double in_exponential_constant)
    : polynomial_coefficients(in_polynomial_coefficients),
      sin_coefficients(in_sin_coefficients),
      cos_coefficients(in_cos_coefficients),
      trig_scale(in_trig_scale),
      reference_density(in_reference_density) {
  if (std::isnan(in_exponential_constant)) {
    // In practice we should always know this at compile time,
    // it should probably be templated on
    has_exponential_prefactor = true;
    exponential_external_constant = in_exponential_constant;
  } else {
    has_exponential_prefactor = false;
    // Maybe this should be a signaling nan?
    exponential_external_constant = std::numeric_limits<double>::quiet_NaN();
  }
}

// Given an expansion h(x) = sum_i f_i(x) , compute int_a^x sum_i f_i(x) e^x +
// F(a)
Enthalpy::Coefficients Enthalpy::Coefficients::compute_exponential_integral(
    std::pair<double, double> initial_condition) {
  // Not yet implemented
  // Doing the exponential integral of something already with an expoenetial
  // prefactor shouldn't be happening
  if (has_exponential_prefactor) {
    ERROR(
        "Attempting to exponentially integrate an EoS decomposition series "
        "with an exponential prefactor!");
  }
  std::vector<double> integral_poly_coeffs(polynomial_coefficients.size() + 1);
  std::vector<double> integral_sin_coeffs(sin_coefficients.size());
  std::vector<double> integral_cos_coeffs(cos_coefficients.size());
  Enthalpy::Coefficients exponential_integral_coefficients = *this;
  // Preprocessing to put coefficients in better basis c_i z^i  = c_i' z^i/i!
  std::vector<double> taylor_series_coefficients(
      polynomial_coefficients.size());
  for (size_t i = 0; i < integral_poly_coeffs.size(); i++) {
    taylor_series_coefficients[i] =
        polynomial_coefficients[i] * static_cast<double>(factorial(i + 1));
  }

  // i indexes elements of the derivative, r sum over lower degree terms
  for (size_t i = 0; i < integral_poly_coeffs.size(); i++) {
    for (size_t r = 0; r <= i; r++) {
      integral_poly_coeffs[i] +=
          ((i - r) % 2 == 0 ? 1 : -1) * taylor_series_coefficients[r];
    }
    // restore the default normalization for the coefficient
    integral_poly_coeffs[i] /= static_cast<double>(factorial(i + 1));
  }
  for (size_t j = 0; j < sin_coefficients.size(); j++) {
    // contribution from the sine terms
    double k = trig_scale;
    integral_sin_coeffs[j] += 1.0 / (square(j) * square(k) + 1.0);
    integral_cos_coeffs[j] -= (j * k) / (square(j) * square(k) + 1.0);
    // contribution from the cosine terms
    integral_cos_coeffs[j] += 1.0 / (square(j) * square(k) + 1.0);
    integral_sin_coeffs[j] += (j * k) / (square(j) * square(k) + 1.0);
  }
  exponential_integral_coefficients.has_exponential_prefactor = true;
  exponential_integral_coefficients.exponential_external_constant = 0.0;
  exponential_integral_coefficients.polynomial_coefficients =
      integral_poly_coeffs;
  exponential_integral_coefficients.sin_coefficients = integral_sin_coeffs;
  exponential_integral_coefficients.cos_coefficients = integral_cos_coeffs;
  double new_constant = std::get<1>(initial_condition) -
                        evaluate_coefficients(exponential_integral_coefficients,
                                              std::get<0>(initial_condition));
  exponential_integral_coefficients.exponential_external_constant =
      new_constant;
  return exponential_integral_coefficients;
}
Enthalpy::Coefficients Enthalpy::Coefficients::compute_derivative() {
  // Slow/Bad?
  std::vector<double> derivative_poly_coeffs(polynomial_coefficients.size());
  std::vector<double> derivative_sin_coeffs(sin_coefficients.size());
  std::vector<double> derivative_cos_coeffs(cos_coefficients.size());
  Enthalpy::Coefficients derivative_coefficients = *this;
  // d/dz e^z \sum_i a_i f(z) = e^z \sum_i a_i f_i(z)  + e^z \sum_i a_i f_i'(z)
  if (has_exponential_prefactor) {
    derivative_poly_coeffs = polynomial_coefficients;
    derivative_sin_coeffs = sin_coefficients;
    derivative_cos_coeffs = cos_coefficients;
    derivative_poly_coeffs[polynomial_coefficients.size() - 1] = 0.0;
    for (size_t i = 0; i < polynomial_coefficients.size() - 1; i++) {
      derivative_poly_coeffs[i] +=
          polynomial_coefficients[i + 1] * static_cast<double>(i + 1);
    }
    for (size_t j = 0; j < sin_coefficients.size(); j++) {
      derivative_cos_coeffs[j] +=
          -sin_coefficients[j] * (static_cast<double>(j + 1) * trig_scale);
      derivative_sin_coeffs[j] +=
          cos_coefficients[j] * (static_cast<double>(j + 1) * trig_scale);
    }
  } else {
    // The final coefficient will be zero because nothing differentiates to it
    derivative_poly_coeffs[polynomial_coefficients.size() - 1] = 0.0;
    for (size_t i = 0; i < polynomial_coefficients.size() - 1; i++) {
      derivative_poly_coeffs[i] =
          polynomial_coefficients[i + 1] * static_cast<double>(i + 1);
    }
    for (size_t j = 0; j < sin_coefficients.size(); j++) {
      derivative_cos_coeffs[j] =
          -sin_coefficients[j] * (static_cast<double>(j + 1) * trig_scale);
      derivative_sin_coeffs[j] =
          cos_coefficients[j] * (static_cast<double>(j + 1) * trig_scale);
    }
  }

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
                   const EquationsOfState::Spectral& lower_spectral)
    : reference_density_(reference_density),
      minimum_density_(min_density),
      maximum_density_(max_density),
      lower_spectral_(lower_spectral),
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
}
double Enthalpy::density_from_x(const double x) const {
  return reference_density_ * exp(x);
}

double Enthalpy::energy_density_from_log_density(const double x) const {
  return reference_density_ *
         evaluate_coefficients(exponential_integral_coefficients_, x);
}

// Evaluate the function represented by the coefficinets at x  = log(rho)
double Enthalpy::evaluate_coefficients(
    const Enthalpy::Coefficients& coefficients, const double x) {
  // Not as good as Horner's rule, so don't use for polynomials
  // The basis_function should be a family of functions indexed by
  // an integer
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
      +[](double z, int n) { return sin((n + 1) * z); }, x,
      coefficients.sin_coefficients, 0.0);
  const double cos_contribution = evaluate_with_function_basis(
      +[](double z, int n) { return cos((n + 1) * z); }, x,
      coefficients.sin_coefficients, 0.0);
  double value =
      polynomial_contribution + (sin_contribution + cos_contribution);
  if (coefficients.has_exponential_prefactor) {
    value *= exp(x);
    value += coefficients.exponential_external_constant;
  }
  return value;
}

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
