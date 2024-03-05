// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp2::points {
/// A line segment extending from `Begin` to `End`, containing `NumberOfPoints`
/// uniformly-spaced points including the endpoints.
template <size_t Dim>
struct LineSegment : public tt::ConformsTo<protocols::Points> {
  using tags_on_target = tmpl::list<>;
  using points_volume_compute_tags = tmpl::list<>;

  struct Begin {
    using type = std::array<double, Dim>;
    static constexpr Options::String help = {"Beginning endpoint"};
  };
  struct End {
    using type = std::array<double, Dim>;
    static constexpr Options::String help = {"Ending endpoint"};
  };
  struct NumberOfPoints {
    using type = size_t;
    static constexpr Options::String help = {
        "Number of points including endpoints"};
    static type lower_bound() { return 2; }
  };
  using options = tmpl::list<Begin, End, NumberOfPoints>;
  static constexpr Options::String help = {
      "A line segment extending from Begin to End, containing "
      "NumberOfPoints uniformly-spaced points including the endpoints."};

  LineSegment(const std::array<double, Dim>& begin,
              const std::array<double, Dim>& end, size_t number_of_points);

  LineSegment() = default;

  /// @{
  /*!
   * \brief Methods specific to the Sphere target that return the input
   * parameters.
   *
   * \note We don't offer the number of points because that can be easily
   * retrieved from the tensor of points.
   */
  const std::array<double, Dim>& begin() const;
  const std::array<double, Dim>& end() const;
  /// @}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
      const;

 private:
  std::array<double, Dim> begin_{};
  std::array<double, Dim> end_{};
  tnsr::I<DataVector, Dim, Frame::NoFrame> points_{};

  template <size_t LocalDim>
  friend bool operator==(const LineSegment<LocalDim>& lhs,
                         const LineSegment<LocalDim>& rhs);
};

template <size_t Dim>
bool operator!=(const LineSegment<Dim>& lhs, const LineSegment<Dim>& rhs);
}  // namespace intrp2::points
