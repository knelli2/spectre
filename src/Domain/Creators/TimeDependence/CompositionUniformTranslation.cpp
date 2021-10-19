// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependence/CompositionUniformTranslation.hpp"

#include <unordered_map>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/FunctionsOfTime/CombineFunctionsOfTime.hpp"
#include "Utilities/CloneUniquePtrs.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace domain::creators::time_dependence {

template <size_t MeshDim>
CompositionUniformTranslation<MeshDim>::CompositionUniformTranslation(
    const std::unique_ptr<TimeDependence<MeshDim>>& uniform_translation0,
    const std::unique_ptr<TimeDependence<MeshDim>>& uniform_translation1)
    : uniform_translation0_(uniform_translation0->get_clone()),
      uniform_translation1_(uniform_translation1->get_clone()),
      coord_map_(domain::push_back(
          dynamic_cast<UniformTranslation<MeshDim>&>(*uniform_translation0_)
              .map_for_composition(),
          dynamic_cast<UniformTranslation<MeshDim, 1>&>(*uniform_translation1_)
              .map_for_composition())) {}

template <size_t MeshDim>
auto CompositionUniformTranslation<MeshDim>::get_clone() const
    -> std::unique_ptr<TimeDependence<MeshDim>> {
  return std::make_unique<CompositionUniformTranslation>(uniform_translation0_,
                                                         uniform_translation1_);
}

template <size_t MeshDim>
auto CompositionUniformTranslation<MeshDim>::block_maps(size_t number_of_blocks)
    const -> std::vector<std::unique_ptr<
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
auto CompositionUniformTranslation<MeshDim>::functions_of_time(
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const -> std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> {
  const auto ut0_f_of_t =
      dynamic_cast<UniformTranslation<MeshDim>&>(*uniform_translation0_)
          .functions_of_time(initial_expiration_times);
  const auto ut1_f_of_t =
      dynamic_cast<UniformTranslation<MeshDim, 1>&>(*uniform_translation1_)
          .functions_of_time(initial_expiration_times);

  return domain::FunctionsOfTime::combine_functions_of_time(ut0_f_of_t,
                                                            ut1_f_of_t);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) \
  template class CompositionUniformTranslation<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef DIM
#undef INSTANTIATION

}  // namespace domain::creators::time_dependence
