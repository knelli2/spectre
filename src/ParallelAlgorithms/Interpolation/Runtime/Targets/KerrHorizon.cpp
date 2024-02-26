// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/Interpolation/Runtime/Targets/KerrHorizon.hpp"

#include <array>
#include <cstddef>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Transpose.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/StrahlkorperFunctions.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrHorizon.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"
#include "Utilities/StdArrayHelpers.hpp"

namespace intrp2::Targets {
KerrHorizon::KerrHorizon(const size_t l_max,
                         const std::array<double, 3>& center, const double mass,
                         const std::array<double, 3>& dimensionless_spin,
                         const intrp::AngularOrdering angular_ordering,
                         const Options::Context& context)
    : l_max_(l_max),
      center_(center),
      mass_(mass),
      dimensionless_spin_(dimensionless_spin),
      angular_ordering_(angular_ordering) {
  if (mass_ <= 0.0) {
    // Check here, rather than put a lower_bound on the Tag, because
    // we want to exclude mass being exactly zero.
    PARSE_ERROR(context, "KerrHorizon expects mass > 0, not " << mass_);
  }
  if (magnitude(dimensionless_spin_) > 1.0) {
    PARSE_ERROR(context, "KerrHorizon expects |dimensionless_spin|<=1, not "
                             << magnitude(dimensionless_spin_));
  }

  // Make a Strahlkorper with the correct shape.
  ylm::Strahlkorper<Frame::NoFrame> strahlkorper{
      l_max_, l_max_,
      get(gr::Solutions::kerr_horizon_radius(
          ::ylm::Spherepack(l_max_, l_max_).theta_phi_points(), mass_,
          dimensionless_spin_)),
      center_};

  ylm::cartesian_coords(make_not_null(&points_), strahlkorper,
                        ylm::radius(strahlkorper),
                        ylm::rhat(ylm::theta_phi(strahlkorper)));

  if (angular_ordering_ == intrp::AngularOrdering::Cce) {
    const auto physical_extents =
        strahlkorper.ylm_spherepack().physical_extents();
    for (size_t i = 0; i < 3; i++) {
      // Can't use not_null version because points get overwritten before they
      // are used
      points_.get(i) =
          transpose(points_.get(i), physical_extents[0], physical_extents[1]);
    }
  }
}

const std::string& KerrHorizon::name() const { return name_; }

size_t KerrHorizon::l_max() const { return l_max_; }
const std::array<double, 3>& KerrHorizon::center() const { return center_; }
double KerrHorizon::mass() const { return mass_; }
const std::array<double, 3>& KerrHorizon::dimensionless_spin() const {
  return dimensionless_spin_;
}
intrp::AngularOrdering KerrHorizon::angular_ordering() const {
  return angular_ordering_;
}

void KerrHorizon::pup(PUP::er& p) {
  Target<3>::pup(p);
  p | name_;
  p | l_max_;
  p | center_;
  p | mass_;
  p | dimensionless_spin_;
  p | angular_ordering_;
  p | points_;
}

const tnsr::I<DataVector, 3, Frame::NoFrame>&
KerrHorizon::target_points_no_frame() const {
  return points_;
}

PUP::able::PUP_ID KerrHorizon::my_PUP_ID = 0;  // NOLINT

bool operator==(const KerrHorizon& lhs, const KerrHorizon& rhs) {
  return lhs.name_ == rhs.name_ and lhs.l_max_ == rhs.l_max_ and
         lhs.center_ == rhs.center_ and lhs.mass_ == rhs.mass_ and
         lhs.dimensionless_spin_ == rhs.dimensionless_spin_ and
         lhs.angular_ordering_ == rhs.angular_ordering_ and
         lhs.points_ == rhs.points_;
}

bool operator!=(const KerrHorizon& lhs, const KerrHorizon& rhs) {
  return not(lhs == rhs);
}
}  // namespace intrp2::Targets
