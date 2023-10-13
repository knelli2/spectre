// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Target.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp::Targets {
/// A line segment extending from `Begin` to `End`, containing `NumberOfPoints`
/// uniformly-spaced points including the endpoints.
template <size_t Dim>
struct LineSegment : Target<Dim> {
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

  /// \cond
  explicit LineSegment(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(Target<Dim>, LineSegment<Dim>);
  /// \endcond

  const std::string& name() const override;

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
  void pup(PUP::er& p) override;

 private:
  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
      const override;

  std::string name_{"LineSegment"};
  std::array<double, Dim> begin_{};
  std::array<double, Dim> end_{};
  tnsr::I<DataVector, Dim, Frame::NoFrame> points_{};

  template <size_t LocalDim>
  friend bool operator==(const LineSegment<LocalDim>& lhs,
                         const LineSegment<LocalDim>& rhs);
};

template <size_t Dim>
bool operator!=(const LineSegment<Dim>& lhs, const LineSegment<Dim>& rhs);
}  // namespace intrp::Targets
