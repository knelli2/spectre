// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"
#include "PointwiseFunctions/Hydro/TagsDeclarations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace hydro {
template <typename DataType, size_t Dim, typename Frame>
void magnetic_pressure(gsl::not_null<Scalar<DataType>*> result,
                       const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
                       const tnsr::I<DataType, Dim, Frame>& magnetic_field,
                       const Scalar<DataType>& lorentz_factor,
                       const tnsr::I<DataType, Dim, Frame>& spatial_velocity);

template <typename DataType, size_t Dim, typename Frame>
Scalar<DataType> magnetic_pressure(
    const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
    const tnsr::I<DataType, Dim, Frame>& magnetic_field,
    const Scalar<DataType>& lorentz_factor,
    const tnsr::I<DataType, Dim, Frame>& spatial_velocity);

namespace Tags {
template <typename DataType, size_t Dim, typename Frame>
struct MagneticPressureCompute : MagneticPressure<DataType>, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::SpatialMetric<Dim, Frame, DataType>,
                 MagneticField<DataType, Dim, Frame>, LorentzFactor<DataType>,
                 SpatialVelocity<DataType, Dim, Frame>>;

  using return_type = Scalar<DataType>;

  static void function(const gsl::not_null<Scalar<DataType>*> result,
                       const tnsr::ii<DataType, Dim, Frame>& spatial_metric,
                       const tnsr::I<DataType, Dim, Frame>& magnetic_field,
                       const Scalar<DataType>& lorentz_factor,
                       const tnsr::I<DataType, Dim, Frame>& spatial_velocity) {
    magnetic_pressure(result, spatial_metric, magnetic_field, lorentz_factor,
                      spatial_velocity);
  }

  using base = MagneticPressure<DataType>;
};
}  // namespace Tags
}  // namespace hydro
