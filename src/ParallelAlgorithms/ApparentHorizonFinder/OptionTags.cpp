// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/ApparentHorizonFinder/OptionTags.hpp"

#include <string>

#include "Domain/Structure/ObjectLabel.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "Options/ParseError.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"
#include "Utilities/Algorithm.hpp"

namespace ah {
HorizonFinderOptions::HorizonFinderOptions(
    HorizonFinderOptions::InitialGuess initial_guess_in,
    const domain::ObjectLabel label, const std::string& frame,
    ::FastFlow fast_flow_in, ::Verbosity verbosity_in,
    const Options::Context& context)
    : object_label(label),
      fast_flow(std::move(fast_flow_in)),
      verbosity(std::move(verbosity_in)) {
  const std::unordered_set<std::string> frames{"Grid", "Distorted", "Inertial"};
  if (not alg::found(frames, frame)) {
    PARSE_ERROR(context, "Frame '"
                             << frame
                             << "' is not a valid frame. Valid frames are "
                             << frames);
  }

  if (frame == "Grid") {
    initial_guess = ylm::Strahlkorper<Frame::Grid>{initial_guess_in.l_max,
                                                   initial_guess_in.radius,
                                                   initial_guess_in.center};
  } else if (frame == "Distorted") {
    initial_guess = ylm::Strahlkorper<Frame::Distorted>{
        initial_guess_in.l_max, initial_guess_in.radius,
        initial_guess_in.center};
  } else {
    initial_guess = ylm::Strahlkorper<Frame::Inertial>{initial_guess_in.l_max,
                                                       initial_guess_in.radius,
                                                       initial_guess_in.center};
  }
}

bool operator==(const HorizonFinderOptions& lhs,
                const HorizonFinderOptions& rhs) {
  return lhs.initial_guess == rhs.initial_guess and
         lhs.object_label == rhs.object_label and
         lhs.fast_flow == rhs.fast_flow and lhs.verbosity == rhs.verbosity;
}
bool operator!=(const HorizonFinderOptions& lhs,
                const HorizonFinderOptions& rhs) {
  return not(lhs == rhs);
}
}  // namespace ah
