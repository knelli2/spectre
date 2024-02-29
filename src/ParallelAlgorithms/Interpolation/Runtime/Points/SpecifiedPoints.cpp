// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Runtime/Points/SpecifiedPoints.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace intrp2::points {

template <size_t Dim>
SpecifiedPoints<Dim>::SpecifiedPoints(
    const std::vector<std::array<double, Dim>>& points) {
  set_number_of_grid_points(make_not_null(&points_), points.size());
  for (size_t d = 0; d < Dim; d++) {
    for (size_t i = 0; i < points.size(); ++i) {
      points_.get(d)[i] = gsl::at(points[i], d);
    }
  }
}

template <size_t Dim>
void SpecifiedPoints<Dim>::pup(PUP::er& p) {
  p | points_;
}

template <size_t Dim>
const tnsr::I<DataVector, Dim, Frame::NoFrame>&
SpecifiedPoints<Dim>::target_points_no_frame() const {
  return points_;
}

template <size_t Dim>
bool operator==(const SpecifiedPoints<Dim>& lhs,
                const SpecifiedPoints<Dim>& rhs) {
  return lhs.points_ == rhs.points_;
}

template <size_t Dim>
bool operator!=(const SpecifiedPoints<Dim>& lhs,
                const SpecifiedPoints<Dim>& rhs) {
  return not(lhs == rhs);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                   \
  template struct SpecifiedPoints<DIM(data)>;                  \
  template bool operator==(const SpecifiedPoints<DIM(data)>&,  \
                           const SpecifiedPoints<DIM(data)>&); \
  template bool operator!=(const SpecifiedPoints<DIM(data)>&,  \
                           const SpecifiedPoints<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef DIM
#undef INSTANTIATE
}  // namespace intrp2::points
