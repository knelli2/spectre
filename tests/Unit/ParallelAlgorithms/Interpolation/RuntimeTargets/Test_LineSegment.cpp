// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/ParallelAlgorithms/Interpolation/RuntimeTargetHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/LineSegment.hpp"

namespace intrp::Targets {
namespace {
template <size_t Dim>
void test_line_segment() {
  const std::string extra_zeros =
      Dim == 1 ? "" : (Dim == 2 ? ", 0.0" : ", 0.0, 0.0");

  const std::string expected_name = "LineSegment";
  const std::string option_string =
      "LineSegment:\n"
      " Begin: [0.0" +
      extra_zeros +
      "]\n"
      " End: [1.0" +
      extra_zeros +
      "]\n"
      " NumberOfPoints: 3\n";

  tnsr::I<DataVector, Dim, Frame::NoFrame> expected_points_no_frame{3_st, 0.0};
  get<0>(expected_points_no_frame)[1] = 0.5;
  get<0>(expected_points_no_frame)[2] = 1.0;

  const auto target = TestHelpers::test_target_all_frames<LineSegment<Dim>>(
      option_string, expected_points_no_frame, expected_name);
  const LineSegment<Dim>& line_segment =
      dynamic_cast<const LineSegment<Dim>&>(*target.get());

  const auto begin = make_array<Dim, double>(0.0);
  auto end = make_array<Dim, double>(0.0);
  end[0] = 1.0;

  CHECK(line_segment.begin() == begin);
  CHECK(line_segment.end() == end);
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.RuntimeTargets.LineSegment",
    "[Unit]") {
  test_line_segment<1>();
  test_line_segment<2>();
  test_line_segment<3>();
}
}  // namespace intrp::Targets
