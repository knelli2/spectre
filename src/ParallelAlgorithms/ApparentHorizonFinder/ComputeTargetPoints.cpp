// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/ApparentHorizonFinder/ComputeTargetPoints.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/StrahlkorperFunctions.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"

namespace ah {
tnsr::I<DataVector, 3, Frame::NoFrame> compute_target_points(
    const ::Strahlkorper<Frame::NoFrame>& old_horizon,
    const ::FastFlow& fast_flow) {
  const size_t L_mesh = fast_flow.current_l_mesh(old_horizon);
  const auto new_horizon =
      ylm::Strahlkorper<Frame::NoFrame>(L_mesh, L_mesh, old_horizon);

  Variables<
      tmpl::list<::Tags::Tempi<0, 2, ::Frame::Spherical<Frame::NoFrame>>,
                 ::Tags::Tempi<1, 3, Frame::NoFrame>, ::Tags::TempScalar<2>>>
      temp_buffer(new_horizon.ylm_spherepack().physical_size());

  auto& theta_phi =
      get<::Tags::Tempi<0, 2, ::Frame::Spherical<Frame::NoFrame>>>(temp_buffer);
  auto& r_hat = get<::Tags::Tempi<1, 3, Frame::NoFrame>>(temp_buffer);
  auto& radius = get<::Tags::TempScalar<2>>(temp_buffer);
  ylm::theta_phi(make_not_null(&theta_phi), new_horizon);
  ylm::rhat(make_not_null(&r_hat), theta_phi);
  ylm::radius(make_not_null(&radius), new_horizon);

  tnsr::I<DataVector, 3, Frame::NoFrame> new_coords{};
  ylm::Tags::CartesianCoordsCompute<Frame::NoFrame>::function(
      make_not_null(&new_coords), new_horizon, radius, r_hat);

  return new_coords;
}
}  // namespace ah
