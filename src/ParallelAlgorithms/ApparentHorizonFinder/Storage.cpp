// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/ApparentHorizonFinder/Storage.hpp"

#include <pup.h>
#include <pup_stl.h>

namespace ah::Storage {
void NumberAndId::pup(PUP::er& p) {
  p | number;
  p | id;
}

bool operator==(const NumberAndId& lhs, const NumberAndId& rhs) {
  return lhs.number == rhs.number;
}
bool operator<(const NumberAndId& lhs, const NumberAndId& rhs) {
  return lhs.number < rhs.number;
}

void VolumeVars::pup(PUP::er& p) {
  p | mesh;
  p | volume_vars_for_horizon_finding;
  p | extra_volume_vars;
}

void VolumeDataAndCallbacks::pup(PUP::er& p) {
  p | frame;
  p | destination;
  p | error_on_failure;
  p | volume_vars_per_element;
  p | callbacks;
}

void InterpolatedVars::pup(PUP::er& p) {
  p | vars;
  p | interpolation_is_done_for_these_elements;
  p | indicies_interpolated_to_thus_far;
}

void SurfaceAndPoints::pup(PUP::er& p) {
  p | strahlkorper;
  p | frame;
  p | cartesian_coords;
  p | block_logical_coords;
}
}  // namespace ah::Storage
