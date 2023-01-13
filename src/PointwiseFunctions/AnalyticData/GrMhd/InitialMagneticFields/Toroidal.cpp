// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/GrMhd/InitialMagneticFields/Toroidal.hpp"

#include <memory>

#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/InitialMagneticFields/InitialMagneticField.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace grmhd::AnalyticData::InitialMagneticFields {

std::unique_ptr<InitialMagneticField> Toroidal::get_clone() const {
  return std::make_unique<Toroidal>(*this);
}

Toroidal::Toroidal(CkMigrateMessage* msg) : InitialMagneticField(msg) {}

void Toroidal::pup(PUP::er& p) {
  InitialMagneticField::pup(p);
  p | pressure_exponent_;
  p | cutoff_pressure_;
  p | vector_potential_amplitude_;
}

// NOLINTNEXTLINE
PUP::able::PUP_ID Toroidal::my_PUP_ID = 0;

Toroidal::Toroidal(const size_t pressure_exponent, const double cutoff_pressure,
                   const double vector_potential_amplitude)
    : pressure_exponent_(pressure_exponent),
      cutoff_pressure_(cutoff_pressure),
      vector_potential_amplitude_(vector_potential_amplitude) {}

template <typename DataType>
tuples::TaggedTuple<hydro::Tags::MagneticField<DataType, 3>>
Toroidal::variables(const tnsr::I<DataType, 3>& coords,
                    const Scalar<DataType>& pressure,
                    const Scalar<DataType>& sqrt_det_spatial_metric,
                    const tnsr::i<DataType, 3>& deriv_pressure) const {
  auto magnetic_field = make_with_value<tnsr::I<DataType, 3>>(coords, 0.0);
  const size_t num_pts = get_size(get(pressure));

  for (size_t i = 0; i < num_pts; ++i) {
    const double pressure_i = get_element(get(pressure), i);
    if (pressure_i < cutoff_pressure_) {
      get_element(magnetic_field.get(0), i) = 0.0;
      get_element(magnetic_field.get(1), i) = 0.0;
      get_element(magnetic_field.get(2), i) = 0.0;
      continue;
    }

    // (p - p_c)^{n_s}
    const double pressure_term =
        2.0 * pow(pressure_i - cutoff_pressure_, pressure_exponent_);
    // n_s * (p - p_c)^{n_s-1}
    const double n_times_pressure_to_n_minus_1 =
        pressure_exponent_ * pow(pressure_i - cutoff_pressure_,
                                 static_cast<int>(pressure_exponent_) - 1);

    const double x = get_element(coords.get(0), i);
    const double y = get_element(coords.get(1), i);
    const double varpi_squared = square(x) + square(y);

    const auto& dp_dx = get_element(deriv_pressure.get(0), i);
    const auto& dp_dy = get_element(deriv_pressure.get(1), i);

    // Assign Bx, By
    get_element(magnetic_field.get(0), i) =
        y * pressure_term +
        varpi_squared * n_times_pressure_to_n_minus_1 * dp_dy;
    get_element(magnetic_field.get(1), i) =
        -(x * pressure_term +
          varpi_squared * n_times_pressure_to_n_minus_1 * dp_dx);
  }

  magnetic_field.get(0) *=
      vector_potential_amplitude_ / get(sqrt_det_spatial_metric);
  magnetic_field.get(1) *=
      vector_potential_amplitude_ / get(sqrt_det_spatial_metric);

  return {std::move(magnetic_field)};
}

bool operator==(const Toroidal& lhs, const Toroidal& rhs) {
  return lhs.pressure_exponent_ == rhs.pressure_exponent_ and
         lhs.cutoff_pressure_ == rhs.cutoff_pressure_ and
         lhs.vector_potential_amplitude_ == rhs.vector_potential_amplitude_;
}

bool operator!=(const Toroidal& lhs, const Toroidal& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                               \
  template tuples::TaggedTuple<hydro::Tags::MagneticField<DTYPE(data), 3>> \
  Toroidal::variables<DTYPE(data)>(                                        \
      const tnsr::I<DTYPE(data), 3>& coords,                               \
      const Scalar<DTYPE(data)>& pressure,                                 \
      const Scalar<DTYPE(data)>& sqrt_det_spatial_metric,                  \
      const tnsr::i<DTYPE(data), 3>& deriv_pressure) const;

GENERATE_INSTANTIATIONS(INSTANTIATE, (double, DataVector))

#undef INSTANTIATE
#undef DTYPE

}  // namespace grmhd::AnalyticData::InitialMagneticFields