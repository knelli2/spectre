// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/Runtime/PointsHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/LineSegment.hpp"
#include "Utilities/MakeArray.hpp"

namespace {
template <size_t Dim>
void test_line_segment() {
  const auto extra_string = [](const std::string& number) -> std::string {
    std::string result = Dim > 1 ? (", " + number) : "";
    if constexpr (Dim > 2) {
      result += ", " + number;
    }
    return result;
  };
  const std::string option_string{"Begin: [1.0" + extra_string("1.0") +
                                  "]\n"
                                  "End: [2.4" +
                                  extra_string("2.4") +
                                  "]\n"
                                  "NumberOfPoints: 15"};

  const size_t n_pts = 15;
  tnsr::I<DataVector, Dim, Frame::NoFrame> expected_points(n_pts);
  for (size_t d = 0; d < Dim; ++d) {
    for (size_t i = 0; i < n_pts; ++i) {
      expected_points.get(d)[i] = 1.0 + 0.1 * i;  // Worked out by hand.
    }
  }

  const auto line_segment =
      RuntimeIntrpTestHelpers::test_points<intrp2::points::LineSegment<Dim>>(
          option_string, expected_points);

  CHECK(line_segment.begin_point() == make_array<Dim, double>(1.0));
  CHECK(line_segment.end_point() == make_array<Dim, double>(2.4));
  CHECK(line_segment.number_of_sets_of_points() == 1);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.Points.LineSegment",
    "[Unit]") {
  test_line_segment<1>();
  test_line_segment<2>();
  test_line_segment<3>();
}
