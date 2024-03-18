// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"

namespace ylm {
/// \ingroup SurfacesGroup
/// \brief Returns a list of `ylm::Strahlkorper`s constructed from reading in
/// spherical harmonic data for a surface at a requested list of times
///
/// \details The `ylm::Strahlkorper`s are constructed by reading in data from
/// an H5 subfile that is expected to be in the format described by
/// `intrp::callbacks::ObserveSurfaceData`. It is assumed that
/// \f$l_{max} = m_{max}\f$.
///
/// \param file_name name of the H5 file containing the surface's spherical
/// harmonic data
/// \param surface_subfile_name name of the subfile (with no leading slash nor
/// the `.dat` extension) within `file_name` that contains the surface's
/// spherical harmonic data to read in
/// \param requested_number_of_times_from_end the number of times to read in
/// starting backwards from the final time found in `surface_subfile_name`
template <typename Frame>
std::vector<ylm::Strahlkorper<Frame>> read_surface_ylm(
    const std::string& file_name, const std::string& surface_subfile_name,
    size_t requested_number_of_times_from_end);

template <typename Frame>
ylm::Strahlkorper<Frame> read_surface_ylm_single_time(
    const std::string& file_name, const std::string& surface_subfile_name,
    double time, double epsilon);
}  // namespace ylm
