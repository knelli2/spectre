// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Ringdown/WrapFillYlmData.hpp"

#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveSurfaceData.hpp"
#include "Utilities/Gsl.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace evolution::Ringdown {
std::vector<double> wrap_fill_ylm_data(const DataVector& coefs,
                                       const double match_time,
                                       const std::array<double, 3>& center,
                                       const size_t l_max) {
  const ylm::Strahlkorper<Frame::Distorted> other{l_max, l_max, 1.0, center};
  const ylm::Strahlkorper<Frame::Distorted> strahlkorper{coefs, other};
  std::vector<std::string> legend{};
  std::vector<double> data{};
  intrp::callbacks::detail::fill_ylm_legend_and_data(
      make_not_null(&legend), make_not_null(&data), strahlkorper, match_time,
      l_max);

  return data;
}
}  // namespace evolution::Ringdown
