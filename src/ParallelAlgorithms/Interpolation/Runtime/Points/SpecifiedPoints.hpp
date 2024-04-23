// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
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
/// A list of specified points to interpolate to.
template <size_t Dim>
struct SpecifiedPoints : public tt::ConformsTo<protocols::Points> {
  using tags_on_target = tmpl::list<>;
  using points_volume_compute_tags = tmpl::list<>;

  struct Points {
    using type = std::vector<std::array<double, Dim>>;
    static constexpr Options::String help = {"Coordinates of each point"};
  };
  using options = tmpl::list<Points>;
  static constexpr Options::String help = {
      "A list of specified points in frame 'Frame'"};

  SpecifiedPoints(const std::vector<std::array<double, Dim>>& points);

  SpecifiedPoints() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
      const;

  /*!
   * \brief Always returns 1 since these are individual points
   */
  size_t number_of_sets_of_points() const;

 private:
  tnsr::I<DataVector, Dim, Frame::NoFrame> points_{};

  template <size_t LocalDim>
  friend bool operator==(const SpecifiedPoints<LocalDim>& lhs,
                         const SpecifiedPoints<LocalDim>& rhs);
};

template <size_t Dim>
bool operator!=(const SpecifiedPoints<Dim>& lhs,
                const SpecifiedPoints<Dim>& rhs);
}  // namespace intrp2::points
