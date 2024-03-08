// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/Runtime/PointsHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/SpecifiedPoints.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"

namespace {
template <size_t Dim, typename Generator>
void test_specified_points(const gsl::not_null<Generator*> gen) {
  std::uniform_int_distribution<size_t> int_dist{1, 15};
  const size_t num_points = int_dist(*gen);
  tnsr::I<DataVector, Dim, Frame::NoFrame> expected_points(num_points);

  std::string option_string{"Points: ["};
  for (size_t i = 0; i < num_points; i++) {
    option_string += "[";
    for (size_t j = 0; j < Dim; j++) {
      const size_t number = int_dist(*gen);

      option_string += std::to_string(number);
      expected_points.get(j)[i] = double(number);

      if (j != Dim - 1) {
        option_string += ",";
      }
    }
    option_string += "]";
    if (i != num_points - 1) {
      option_string += ",";
    }
  }
  option_string += "]";

  // No need to use the return value
  RuntimeIntrpTestHelpers::test_points<intrp2::points::SpecifiedPoints<Dim>>(
      option_string, expected_points);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.Points.SpecifiedPoints",
    "[Unit]") {
  MAKE_GENERATOR(gen);
  test_specified_points<1>(make_not_null(&gen));
  test_specified_points<2>(make_not_null(&gen));
  test_specified_points<3>(make_not_null(&gen));
}
