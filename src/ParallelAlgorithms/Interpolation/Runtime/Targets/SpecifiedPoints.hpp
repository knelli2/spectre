// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Targets/Target.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp2::Targets {
/// A list of specified points to interpolate to.
template <size_t Dim>
struct SpecifiedPoints : Target<Dim> {
  struct Points {
    using type = std::vector<std::array<double, Dim>>;
    static constexpr Options::String help = {"Coordinates of each point"};
  };
  using options = tmpl::list<Points>;
  static constexpr Options::String help = {
      "A list of specified points in frame 'Frame'"};

  SpecifiedPoints(const std::vector<std::array<double, Dim>>& points);

  SpecifiedPoints() = default;

  /// \cond
  explicit SpecifiedPoints(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(Target<Dim>, SpecifiedPoints<Dim>);
  /// \endcond

  const std::string& name() const override;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
      const override;

  std::string name_{"SpecifiedPoints"};
  tnsr::I<DataVector, Dim, Frame::NoFrame> points_{};

  template <size_t LocalDim>
  friend bool operator==(const SpecifiedPoints<LocalDim>& lhs,
                         const SpecifiedPoints<LocalDim>& rhs);
};

template <size_t Dim>
bool operator!=(const SpecifiedPoints<Dim>& lhs,
                const SpecifiedPoints<Dim>& rhs);
}  // namespace intrp2::Targets
