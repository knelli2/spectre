// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependence/CompositionCubicScaleAndUniformRotationAboutZAxis.hpp"

#include <unordered_map>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/FunctionsOfTime/CombineFunctionsOfTime.hpp"
#include "Utilities/CloneUniquePtrs.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace domain::creators::time_dependence {

template <size_t MeshDim>
CompositionCubicScaleAndUniformRotationAboutZAxis<MeshDim>::
    CompositionCubicScaleAndUniformRotationAboutZAxis(
        const std::unique_ptr<TimeDependence<MeshDim>>& cubic_scale,
        const std::unique_ptr<TimeDependence<MeshDim>>&
            uniform_rotation_about_z_axis)
    : cubic_scale_(cubic_scale->get_clone()),
      uniform_rotation_about_z_axis_(
          uniform_rotation_about_z_axis->get_clone()),
      coord_map_(
          domain::push_back(dynamic_cast<CubicScale<MeshDim>&>(*cubic_scale_)
                                .map_for_composition(),
                            dynamic_cast<UniformRotationAboutZAxis<MeshDim>&>(
                                *uniform_rotation_about_z_axis_)
                                .map_for_composition())) {}

template <size_t MeshDim>
auto CompositionCubicScaleAndUniformRotationAboutZAxis<MeshDim>::get_clone()
    const -> std::unique_ptr<TimeDependence<MeshDim>> {
  return std::make_unique<CompositionCubicScaleAndUniformRotationAboutZAxis>(
      cubic_scale_, uniform_rotation_about_z_axis_);
}

template <size_t MeshDim>
auto CompositionCubicScaleAndUniformRotationAboutZAxis<MeshDim>::block_maps(
    size_t number_of_blocks) const
    -> std::vector<std::unique_ptr<
        domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, MeshDim>>> {
  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, MeshDim>>>
      result{number_of_blocks};
  result[0] = std::make_unique<CoordMap>(coord_map_);
  for (size_t i = 1; i < number_of_blocks; ++i) {
    result[i] = result[0]->get_clone();
  }
  return result;
}

template <size_t MeshDim>
auto CompositionCubicScaleAndUniformRotationAboutZAxis<
    MeshDim>::functions_of_time(const std::unordered_map<std::string, double>&
                                    initial_expiration_times) const
    -> std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> {
  const auto cubic_scale_f_of_t =
      dynamic_cast<CubicScale<MeshDim>&>(*cubic_scale_)
          .functions_of_time(initial_expiration_times);
  const auto rotation_f_of_t =
      dynamic_cast<UniformRotationAboutZAxis<MeshDim>&>(
          *uniform_rotation_about_z_axis_)
          .functions_of_time(initial_expiration_times);

  return domain::FunctionsOfTime::combine_functions_of_time(cubic_scale_f_of_t,
                                                            rotation_f_of_t);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) \
  template class CompositionCubicScaleAndUniformRotationAboutZAxis<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (2, 3))

#undef DIM
#undef INSTANTIATION

}  // namespace domain::creators::time_dependence
