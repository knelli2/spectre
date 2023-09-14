// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <vector>

#include "Domain/Structure/ObjectLabel.hpp"
#include "Framework/TestCreation.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/OptionTags.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace Frame {
struct Grid;
struct Distorted;
struct Inertial;
}  // namespace Frame

namespace ah {
namespace {

std::string make_option_string(domain::ObjectLabel label,
                               const std::string& frame) {
  return "InitialGuess:\n"
         "  LMax: 4\n"
         "  Radius: 2.1\n"
         "  Center: [0.1, 0.2, 0.3]\n"
         "Name: " +
         (label == domain::ObjectLabel::None ? "None" : get_output(label)) +
         "\n"
         "Frame: " +
         frame +
         "\n"
         "FastFlow:\n"
         "  Flow: Fast\n"
         "  Alpha: 1.0\n"
         "  Beta: 0.5\n"
         "  AbsTol: 1e-12\n"
         "  TruncationTol: 1e-2\n"
         "  DivergenceTol: 1.2\n"
         "  DivergenceIter: 5\n"
         "  MaxIts: 100\n"
         "Verbosity: Silent\n";
}

template <typename Frame>
void test_single_ah_option(domain::ObjectLabel label, const Frame& /*frame*/) {
  const HorizonFinderOptions horizons_from_options =
      TestHelpers::test_creation<HorizonFinderOptions>(
          make_option_string(label, pretty_type::name<Frame>()));

  const HorizonFinderOptions expected_options{
      HorizonFinderOptions::InitialGuess{4_st, 2.1, std::array{0.1, 0.2, 0.3}},
      label, pretty_type::name<Frame>(),
      FastFlow{FastFlow::FlowType::Fast, 1.0, 0.5, 1.e-12, 1.e-2, 1.2, 5,
               100_st},
      ::Verbosity::Silent};

  CHECK(expected_options == horizons_from_options);
  CHECK(std::holds_alternative<ylm::Strahlkorper<Frame>>(
      horizons_from_options.initial_guess));
}

void test_ah_options() {
  tmpl::for_each<tmpl::list<Frame::Grid, Frame::Distorted, Frame::Inertial>>(
      [](auto frame_v) {
        using Frame = tmpl::type_from<decltype(frame_v)>;
        for (const domain::ObjectLabel label :
             {domain::ObjectLabel::A, domain::ObjectLabel::B,
              domain::ObjectLabel::C, domain::ObjectLabel::None}) {
          test_single_ah_option(label, Frame{});
        }
      });

  CHECK_THROWS_WITH(
      (TestHelpers::test_creation<HorizonFinderOptions>(
          make_option_string(domain::ObjectLabel::A, "NotAFrame"))),
      Catch::Matchers::ContainsSubstring(
          "Frame 'NotAFrame' is not a valid frame. Valid frames are"));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ApparentHorizonFinder.OptionTags",
                  "[ApparentHorizonFinder][Unit]") {
  test_ah_options();
}
}  // namespace ah
