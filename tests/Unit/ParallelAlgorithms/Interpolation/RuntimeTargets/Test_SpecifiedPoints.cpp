// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/RuntimeTargetHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/SpecifiedPoints.hpp"
#include "Utilities/Gsl.hpp"

namespace intrp::Targets {
namespace {
template <size_t Dim, typename Generator>
void test_specified_points(const gsl::not_null<Generator*> generator) {
  const auto to_string = [](const double d) -> std::string {
    std::stringstream ss{};
    ss << std::setprecision(std::numeric_limits<double>::digits10 + 4)
       << std::scientific << d;
    return ss.str();
  };

  std::uniform_real_distribution<double> value_dist{-1.0, 1.0};
  std::uniform_int_distribution<size_t> size_dist{1, 10};

  const size_t size = size_dist(*generator);
  tnsr::I<DataVector, Dim, Frame::NoFrame> expected_points_no_frame{};
  for (size_t d = 0; d < Dim; d++) {
    expected_points_no_frame.get(d) = make_with_random_values<DataVector>(
        generator, make_not_null(&value_dist), size);
  }

  std::stringstream ss{};
  ss << "[";
  for (size_t i = 0; i < size; i++) {
    if (i != 0) {
      ss << ",";
    }
    ss << "[";
    for (size_t d = 0; d < Dim; d++) {
      if (d != 0) {
        ss << ",";
      }
      ss << to_string(expected_points_no_frame.get(d)[i]);
    }
    ss << "]";
  }
  ss << "]";
  const std::string expected_name = "SpecifiedPoints";
  const std::string option_string =
      "SpecifiedPoints:\n"
      " Points: " +
      ss.str() + "\n";

  TestHelpers::test_target_all_frames<SpecifiedPoints<Dim>>(
      option_string, expected_points_no_frame, expected_name);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.RuntimeTargets.SpecifiedPoints",
    "[Unit]") {
  MAKE_GENERATOR(generator);
  test_specified_points<1>(make_not_null(&generator));
  test_specified_points<2>(make_not_null(&generator));
  test_specified_points<3>(make_not_null(&generator));
}
}  // namespace intrp::Targets
