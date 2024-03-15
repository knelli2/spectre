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
std::vector<DataVector> strahlkorper_coefs_in_ringdown_distorted_frame(
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name, std::vector<double> ahc_times,
    const size_t requested_number_of_times_from_end, const double match_time,
    const double settling_timescale,
    const std::array<double, 3> exp_func_and_2_derivs,
    const std::array<double, 3> exp_outer_bdry_func_and_2_derivs,
    const std::array<std::array<double, 4>, 3> rot_func_and_2_derivs);
}  // namespace evolution::Ringdown
