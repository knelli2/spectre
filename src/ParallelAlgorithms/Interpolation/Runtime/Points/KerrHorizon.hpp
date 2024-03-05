// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/AngularOrdering.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct NoFrame;
}  // namespace Frame
/// \endcond

namespace intrp2::points {
/*!
 * \brief A surface that conforms to the horizon of a Kerr black hole in
 * Kerr-Schild coordinates.
 *
 * \details In addition to the parameters for the Kerr black hole, this holder
 * contains the `LMax` which encodes the angular resolution of the spherical
 * harmonic basis and `AngularOrdering` which encodes the collocation
 * ordering.
 */
struct KerrHorizon : public tt::ConformsTo<protocols::Points> {
  using tags_on_target = tmpl::list<>;
  using points_volume_compute_tags = tmpl::list<>;

  struct LMax {
    using type = size_t;
    static constexpr Options::String help = {
        "KerrHorizon is expanded in Ylms up to l=LMax"};
  };
  struct Center {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {"Center of black hole"};
  };
  struct Mass {
    using type = double;
    static constexpr Options::String help = {"Mass of black hole"};
  };
  struct DimensionlessSpin {
    using type = std::array<double, 3>;
    static constexpr Options::String help = {
        "Dimensionless spin of black hole"};
  };
  struct AngularOrdering {
    using type = intrp2::AngularOrdering;
    static constexpr Options::String help = {
        "Chooses theta,phi ordering in 2d array"};
  };
  using options =
      tmpl::list<LMax, Center, Mass, DimensionlessSpin, AngularOrdering>;
  static constexpr Options::String help = {
      "A Strahlkorper conforming to the horizon (in Kerr-Schild coordinates)"
      " of a Kerr black hole with a specified center, mass, and spin."};

  KerrHorizon(size_t l_max, const std::array<double, 3>& center, double mass,
              const std::array<double, 3>& dimensionless_spin,
              intrp2::AngularOrdering angular_ordering,
              const Options::Context& context = {});

  KerrHorizon() = default;

  /// @{
  /*!
   * \brief Methods specific to the KerrHorizon target that return the input
   * parameters.
   */
  size_t l_max() const;
  const std::array<double, 3>& center() const;
  double mass() const;
  const std::array<double, 3>& dimensionless_spin() const;
  intrp2::AngularOrdering angular_ordering() const;
  /// @}

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  const tnsr::I<DataVector, 3, Frame::NoFrame>& target_points_no_frame() const;

 private:
  size_t l_max_{};
  std::array<double, 3> center_{};
  double mass_{};
  std::array<double, 3> dimensionless_spin_{};
  intrp2::AngularOrdering angular_ordering_{};
  tnsr::I<DataVector, 3, Frame::NoFrame> points_{};

  friend bool operator==(const KerrHorizon& lhs, const KerrHorizon& rhs);
};

bool operator!=(const KerrHorizon& lhs, const KerrHorizon& rhs);
}  // namespace intrp2::points
