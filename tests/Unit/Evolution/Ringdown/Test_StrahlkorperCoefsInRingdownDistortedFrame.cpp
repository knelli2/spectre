// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>

#include "Evolution/Ringdown/StrahlkorperCoefsInRingdownDistortedFrame.hpp"

SPECTRE_TEST_CASE("Unit.Evolution.Ringdown.SCoefsRDis", "[Unit][Evolution]") {
  const std::string& path_to_horizons_h5{
      "/Users/geoffrey/Downloads/spectre_vis_of_the_day/010224/"
      "l10lessvolendearlier/GhBinaryBlackHoleReductionData.h5"};
  const std::string surface_subfile_name{"ObservationAhC_Ylm"};
  const size_t requested_number_of_times_from_end{81};
  const double match_time{4899.799999995201};
  const double settling_timescale{10.0};
  const std::array<double, 3> exp_outer_bdry_func_and_2_derivs{
      {0.9951007101717868, -1.000104099257707e-06, 4.247795367633651e-14}};
  const std::array<std::array<double, 4>, 3> rot_func_and_2_derivs{
      {{0.29654264172111816, 1.1396897376847262e-05, 2.070559307581652e-06,
        0.9550196131529888},
       {-0.13041764054305988, -3.3257819465582044e-07, -1.3652951629165736e-06,
        0.040495913515962566},
       {-0.02082157558287204, -1.3470477437102956e-07, -2.3496973780978276e-07,
        -0.013061715986062053}}};
  evolution::Ringdown::strahlkorper_coefs_in_ringdown_distorted_frame(
      path_to_horizons_h5, surface_subfile_name,
      requested_number_of_times_from_end, match_time, settling_timescale,
      exp_outer_bdry_func_and_2_derivs, rot_func_and_2_derivs);
}
