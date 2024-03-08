// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Framework/TestingFramework.hpp"

#include <concepts>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"

namespace RuntimeIntrpTestHelpers {
template <typename PointsType, size_t Dim>
requires std::equality_comparable<PointsType> PointsType
test_points(const std::string& option_string,
            const tnsr::I<DataVector, Dim, Frame::NoFrame>& expected_points) {
  auto created_points_class =
      TestHelpers::test_creation<PointsType>(option_string);
  test_serialization(created_points_class);

  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points =
      created_points_class.target_points_no_frame();

  CHECK_ITERABLE_APPROX(expected_points, target_points);

  return created_points_class;
}
}  // namespace RuntimeIntrpTestHelpers
