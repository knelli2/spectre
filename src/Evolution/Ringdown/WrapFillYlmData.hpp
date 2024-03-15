// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataVector.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

/// Documentation here
namespace evolution::Ringdown {
std::vector<double> wrap_fill_ylm_data(const DataVector& coefs,
                                       const double match_time,
                                       const std::array<double, 3>& center,
                                       const size_t l_max);
}  // namespace evolution::Ringdown
