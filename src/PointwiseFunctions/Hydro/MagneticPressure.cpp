// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/Hydro/MagneticPressure.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace hydro {

template <typename DataType, size_t Dim, typename Frame>
void magnetic_pressure(gsl::not_null<Scalar<DataType>*> result,
                       const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
                       const tnsr::I<DataType, Dim, Frame>& magnetic_field,
                       const Scalar<DataType>& lorentz_factor,
                       const tnsr::I<DataType, Dim, Frame>& spatial_velocity) {
  destructive_resize_components(result, get_size(get(lorentz_factor)));
  Scalar<DataType> magnetic_field_squared;
  Scalar<DataType> magnetic_field_dot_velocity;
  destructive_resize_components(make_not_null(&magnetic_field_squared),
                                get_size(get(lorentz_factor)));
  destructive_resize_components(make_not_null(&magnetic_field_dot_velocity),
                                get_size(get(lorentz_factor)));

  get(magnetic_field_squared) = 0.0;
  get(magnetic_field_dot_velocity) = 0.0;
  for (size_t i = 0; i < Dim; i++) {
    for (size_t j = 0; j < Dim; j++) {
      get(magnetic_field_squared) += magnetic_field.get(i) *
                                     magnetic_field.get(j) *
                                     spatial_metric.get(i, j);
      get(magnetic_field_dot_velocity) += magnetic_field.get(i) *
                                          spatial_velocity.get(j) *
                                          spatial_metric.get(i, j);
    }
  }

  const DataType square_lorentz = square(get(lorentz_factor));
  const DataType ratio_squares = get(magnetic_field_squared) / square_lorentz;
  const DataType square_b_dot_v = square(get(magnetic_field_dot_velocity));

  // b^2 = B^2 / W^2 + (B^iv_i)^2
  get(*result) = 0.5 * (ratio_squares + square_b_dot_v);
}

template <typename DataType, size_t Dim, typename Frame>
Scalar<DataType> magnetic_pressure(
    const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
    const tnsr::I<DataType, Dim, Frame>& magnetic_field,
    const Scalar<DataType>& lorentz_factor,
    const tnsr::I<DataType, Dim, Frame>& spatial_velocity) {
  Scalar<DataType> result{};
  magnetic_pressure(make_not_null(&result), spatial_metric, magnetic_field,
                    lorentz_factor, spatial_velocity);
  return result;
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DIM(data) BOOST_PP_TUPLE_ELEM(1, data)
#define FR(data) BOOST_PP_TUPLE_ELEM(2, data)

#define INSTANTIATE(_, data)                                              \
  template void magnetic_pressure(                                        \
      const gsl::not_null<Scalar<DTYPE(data)>*> result,                   \
      const tnsr::ii<DTYPE(data), DIM(data), FR(data)>& spatial_metric,   \
      const tnsr::I<DTYPE(data), DIM(data), FR(data)>& magnetic_field,    \
      const Scalar<DTYPE(data)>& lorentz_factor,                          \
      const tnsr::I<DTYPE(data), DIM(data), FR(data)>& spatial_velocity); \
  template Scalar<DTYPE(data)> magnetic_pressure(                         \
      const tnsr::ii<DTYPE(data), DIM(data), FR(data)>& spatial_metric,   \
      const tnsr::I<DTYPE(data), DIM(data), FR(data)>& magnetic_field,    \
      const Scalar<DTYPE(data)>& lorentz_factor,                          \
      const tnsr::I<DTYPE(data), DIM(data), FR(data)>& spatial_velocity);

GENERATE_INSTANTIATIONS(INSTANTIATE, (double, DataVector), (1, 2, 3),
                        (Frame::Grid, Frame::Inertial))

#undef DTYPE
#undef DIM
#undef FR
#undef INSTANTIATE
}  // namespace hydro
