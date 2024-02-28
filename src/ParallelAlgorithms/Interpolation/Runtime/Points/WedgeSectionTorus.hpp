// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp2::points {
/*!
 * \brief A solid torus of points, useful, e.g., when measuring data from an
 * accretion disc.
 *
 * The torus's cross section (e.g., a cut at $\phi=0$) is a wedge-like shape
 * bounded by $r_{\text{min}} \le r \le r_{\text{max}}$ and
 * $\theta_{\text{min}} \le \theta \le \theta_{\text{max}}$.
 *
 * The grid points are located on surfaces of constant $r$, $\theta$, and
 * $\phi$. There are `NumberRadialPoints` points in the radial direction between
 * `MinRadius` and `MaxRadius` (including these endpoints); `NumberThetaPoints`
 * points in the $\theta$ direction between `MinTheta` and `MaxTheta` (including
 * these endpoints); `NumberPhiPoints` points in the $\phi$ direction (with one
 * point always at $\phi=0$).
 *
 * By default, the points follow a Legendre Gauss-Lobatto distribution in the
 * $r$ and $\theta$ directions, and a uniform distribution in the $\phi$
 * direction. The distribution in the $r$ (and/or $\theta$) direction can be
 * made uniform using the `UniformRadialGrid` (and/or `UniformThetaGrid`)
 * option.
 *
 * The `target_points` form a 3D mesh ordered with $r$ varying fastest, then
 * $\theta$, and finally $\phi$ varying slowest.
 */
struct WedgeSectionTorus : protocols::Points {
  using tags_on_target = tmpl::list<>;
  using points_volume_compute_tags = tmpl::list<>;

  struct RadialBounds {
    using type = std::array<double, 2>;
    static constexpr Options::String help = {
        "Inner and outer radius of torus, respectively."};
  };
  struct ThetaBounds {
    using type = std::array<double, 2>;
    static constexpr Options::String help = {
        "Angle of top and bottom of wedge (in radians), respectively."};
  };
  struct NumberOfGridPoints {
    using type = std::array<size_t, 3>;
    static constexpr Options::String help = {
        "Number of grid points in the phi, theta, radial directions, "
        "respectively."};
  };
  struct UniformRadialGrid {
    using type = bool;
    static constexpr Options::String help = {"Use uniform radial grid"};
  };
  struct UniformThetaGrid {
    using type = bool;
    static constexpr Options::String help = {"Use uniform theta grid"};
  };
  using options = tmpl::list<RadialBounds, ThetaBounds, NumberOfGridPoints,
                             UniformRadialGrid, UniformThetaGrid>;
  static constexpr Options::String help = {
      "A torus extending from MinRadius to MaxRadius in r, MinTheta to MaxTheta"
      " in theta, and 2pi in phi."};

  WedgeSectionTorus(const std::array<double, 2>& radial_bounds,
                    const std::array<double, 2>& theta_bounds,
                    const std::array<size_t, 3>& number_of_grid_points,
                    bool use_uniform_radial_grid, bool use_uniform_theta_grid,
                    const Options::Context& context = {});

  WedgeSectionTorus() = default;

  /// @{
  /*!
   * \brief Methods specific to the KerrHorizon target that return the input
   * parameters.
   */
  const std::array<double, 2>& radial_bounds() const;
  const std::array<double, 2>& theta_bounds() const;
  /// The order of the points is phi, theta, radial. Not x, y, z.
  const std::array<size_t, 3>& number_of_grid_points() const;
  bool use_uniform_radial_grid() const;
  bool use_uniform_theta_grid() const;
  /// @}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  const tnsr::I<DataVector, 3, Frame::NoFrame>& target_points_no_frame() const;

 private:
  std::string name_{"WedgeSectionTorus"};
  std::array<double, 2> radial_bounds_{};
  std::array<double, 2> theta_bounds_{};
  std::array<size_t, 3> number_of_grid_points_{};
  bool use_uniform_radial_grid_{};
  bool use_uniform_theta_grid_{};
  tnsr::I<DataVector, 3, Frame::NoFrame> points_{};

  friend bool operator==(const WedgeSectionTorus& lhs,
                         const WedgeSectionTorus& rhs);
};

bool operator!=(const WedgeSectionTorus& lhs, const WedgeSectionTorus& rhs);
}  // namespace intrp2::points
