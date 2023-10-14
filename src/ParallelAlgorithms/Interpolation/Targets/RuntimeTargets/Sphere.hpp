// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <variant>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Target.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp::Targets {
/*!
 * \brief A series of concentric spherical surfaces.
 *
 * \details The parameter `LMax` sets the number of collocation points on
 * each spherical surface equal to `(l_max + 1) * (2 * l_max + 1)`. The
 * parameter `AngularOrdering` encodes the collocation ordering. For example,
 * the apparent horizon finder relies on spherepack routines that require
 * `Strahlkorper` for `AngularOrdering`, and using these surfaces for a CCE
 * worldtube requires `Cce` for `AngularOrdering`.
 */
struct Sphere : Target<3> {
  struct LMax {
    using type = size_t;
    static constexpr Options::String help = {
        "The number of collocation points on each sphere will be equal to "
        "`(l_max + 1) * (2 * l_max + 1)`"};
  };
  struct Center {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {"Center of every sphere"};
  };
  struct Radius {
    using type = std::variant<double, std::vector<double>>;
    static constexpr Options::String help = {"Radius of the sphere(s)"};
  };
  struct AngularOrdering {
    using type = intrp::AngularOrdering;
    static constexpr Options::String help = {
        "Chooses theta,phi ordering in 2d array"};
  };
  using options = tmpl::list<LMax, Center, Radius, AngularOrdering>;
  static constexpr Options::String help = {
      "An arbitrary number of spherical surfaces."};

  Sphere(size_t l_max, const std::array<double, 3>& center,
         const typename Radius::type& input_radii,
         intrp::AngularOrdering angular_ordering,
         const Options::Context& context = {});

  Sphere() = default;

  /// \cond
  explicit Sphere(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(Target, Sphere);
  /// \endcond

  const std::string& name() const override;

  /// @{
  /*!
   * \brief Methods specific to the Sphere target that return the input
   * parameters.
   */
  size_t l_max() const;
  const std::array<double, 3>& center() const;
  const std::set<double>& radii() const;
  intrp::AngularOrdering angular_ordering() const;
  /// @}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  const tnsr::I<DataVector, 3, Frame::NoFrame>& target_points_no_frame()
      const override;

  std::string name_{"Sphere"};
  size_t l_max_{};
  std::array<double, 3> center_{};
  std::set<double> radii_{};
  intrp::AngularOrdering angular_ordering_{};
  tnsr::I<DataVector, 3, Frame::NoFrame> points_{};

  friend bool operator==(const Sphere& lhs, const Sphere& rhs);
};

template <size_t Dim>
bool operator!=(const Sphere& lhs, const Sphere& rhs);
}  // namespace intrp::Targets
