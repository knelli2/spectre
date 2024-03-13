// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/Sphere.hpp"
#include "Domain/Creators/SphereTimeDependentMaps.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/IO/ReadSurfaceYlm.hpp"
#include "Parallel/Printf.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace evolution::Ringdown {
void strahlkorper_coefs_in_ringdown_distorted_frame(
    const std::string& path_to_horizons_h5,
    const std::string& surface_subfile_name,
    const size_t requested_number_of_times_from_end, const double match_time,
    const double settling_timescale,
    const std::array<double, 3> exp_outer_bdry_func_and_2_derivs,
    const std::array<std::array<double, 4>, 3> rot_func_and_2_derivs) {
  const std::vector<ylm::Strahlkorper<Frame::Inertial>>& ahc_inertial_h5 =
      ylm::read_surface_ylm<Frame::Inertial>(
          path_to_horizons_h5, surface_subfile_name,
          requested_number_of_times_from_end);

  const double initial_time{0.0};
  const size_t l_max{ahc_inertial_h5[0].l_max()};
  Parallel::printf("Match_time: %1.15f\n", match_time);

  const domain::creators::sphere::TimeDependentMapOptions::TranslationMapOptions
      translation_map_options{
          {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}}};
  const domain::creators::sphere::TimeDependentMapOptions::ShapeMapOptions
      shape_map_options{l_max, std::nullopt};
  const domain::creators::sphere::TimeDependentMapOptions::ExpansionMapOptions
      expansion_map_options{{0.0, 0.0, 0.0},
                            settling_timescale,
                            exp_outer_bdry_func_and_2_derivs,
                            settling_timescale};
  const domain::creators::sphere::TimeDependentMapOptions::RotationMapOptions
      rotation_map_options{rot_func_and_2_derivs, settling_timescale};
  const domain::creators::sphere::TimeDependentMapOptions
      time_dependent_map_options{initial_time, shape_map_options,
                                 rotation_map_options, expansion_map_options,
                                 translation_map_options};
  Parallel::printf("TEST\n");
//   const domain::creators::Sphere domain_creator{
//       0.01,
//       200.0,
//       domain::creators::Sphere::Excision{nullptr},
//       0,
//       5,
//       std::nullopt,
//       {100.0},
//       domain::CoordinateMaps::Distribution::Linear,
//       ShellWedges::All,
//       time_dependent_map_options};
}
}  // namespace evolution::Ringdown
