// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/LineSegment.hpp"

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace intrp::Targets {

template <size_t Dim>
LineSegment<Dim>::LineSegment(const std::array<double, Dim>& begin,
                              const std::array<double, Dim>& end,
                              const size_t number_of_points)
    : begin_(begin), end_(end) {
  const double fractional_distance =
      1.0 / (static_cast<double>(number_of_points) - 1);
  set_number_of_grid_points(make_not_null(&points_), number_of_points);
  for (size_t n = 0; n < number_of_points; ++n) {
    for (size_t d = 0; d < Dim; ++d) {
      points_.get(d)[n] =
          gsl::at(begin_, d) + static_cast<double>(n) * fractional_distance *
                                   (gsl::at(end_, d) - gsl::at(begin_, d));
    }
  }
}

template <size_t Dim>
const std::string& LineSegment<Dim>::name() const {
  return name_;
}

template <size_t Dim>
const std::array<double, Dim>& LineSegment<Dim>::begin() const {
  return begin_;
}
template <size_t Dim>
const std::array<double, Dim>& LineSegment<Dim>::end() const {
  return end_;
}

template <size_t Dim>
void LineSegment<Dim>::pup(PUP::er& p) {
  Target<Dim>::pup(p);
  p | name_;
  p | begin_;
  p | end_;
  p | points_;
}

template <size_t Dim>
const tnsr::I<DataVector, Dim, Frame::NoFrame>&
LineSegment<Dim>::target_points_no_frame() const {
  return points_;
}

template <size_t Dim>
PUP::able::PUP_ID LineSegment<Dim>::my_PUP_ID = 0;  // NOLINT

template <size_t Dim>
bool operator==(const LineSegment<Dim>& lhs, const LineSegment<Dim>& rhs) {
  return lhs.name_ == rhs.name_ and lhs.begin_ == rhs.begin_ and
         lhs.end_ == rhs.end_ and lhs.points_ == rhs.points_;
}

template <size_t Dim>
bool operator!=(const LineSegment<Dim>& lhs, const LineSegment<Dim>& rhs) {
  return not(lhs == rhs);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                               \
  template struct LineSegment<DIM(data)>;                  \
  template bool operator==(const LineSegment<DIM(data)>&,  \
                           const LineSegment<DIM(data)>&); \
  template bool operator!=(const LineSegment<DIM(data)>&,  \
                           const LineSegment<DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef DIM
#undef INSTANTIATE
}  // namespace intrp::Targets
