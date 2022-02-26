// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cmath>
#include <cstddef>
#include <pup.h>

#include "ApparentHorizons/ObjectLabel.hpp"
#include "ControlSystem/ApparentHorizons/Measurements.hpp"
#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataVector.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackIterator.hpp"
#include "Options/Options.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain::Tags {
template <size_t Dim>
struct Domain;
struct FunctionsOfTime;
}  // namespace domain::Tags
namespace Frame {
struct Distorted;
struct Grid;
}  // namespace Frame
/// \endcond

namespace control_system {
namespace ControlErrors {
namespace detail {
template <::ah::ObjectLabel Horizon>
std::string excision_sphere_name() {
  return Horizon == ::ah::ObjectLabel::A ? "ObjectAExcisionSphere"
                                         : "ObjectBExcisionSphere";
}

template <::ah::ObjectLabel Horizon>
std::string size_name() {
  return Horizon == ::ah::ObjectLabel::A ? "SizeA" : "SizeB";
}

// This is the Y00 shperical harmonic expressed in terms of SPHEREPACK
// coefficients which has an extra factor of sqrt(2/pi) in front of the
// original coefficient sqrt(1/4pi)
SPECTRE_ALWAYS_INLINE double y00_coef() { return sqrt(0.5) / M_PI; }
}  // namespace detail

/*!
 * \brief Control error in the
 * \link domain::CoordinateMaps::TimeDependent::Shape Shape \endlink coordinate
 * map
 *
 * \details Computes the error for each \f$ l,m \f$ mode of the shape map except
 * for \f$ l=0 \f$ and \f$ l=1 \f$. This is because the \f$ l=0 \f$ mode is
 * controlled by charactestic speed control, and the \f$ l=1 \f$ mode is
 * controlled by the \link
 * domain::CoordinateMaps::TimeDependent::Translation Translation \endlink map.
 * The equation for the control error is Eq. (77) from \cite Hemberger2012jz
 * which is
 *
 * \f[ Q_{lm} = -\frac{r_\mathrm{EB} - Y_{00}\lambda_{00}(t)}{Y_{00}S_{00}}
 * S_{lm} - \lambda_{lm}(t), \quad l>1
 * \f]
 *
 * where \f$ r_\mathrm{EB} \f$ is the radius of the excision boundary in the
 * grid frame, \f$ \lambda_{00}(t) \f$ is the size map parameter, and \f$ S_{lm}
 * \f$ are the coefficients of the harmonic expansion of the apparent horizon.
 *
 * Requirements:
 * - This control error requires that there be at least one excision surface in
 *   the simulation
 * - Currently this control system can only be used with the \link
 *   control_system::Systems::Shape Shape \endlink control system
 *
 * \note The map parameters \f$ \lambda_{lm}(t) \f$ and the coefficients of the
 * apparent horizon \f$ S_{lm} \f$ are stored as SPHEREPACK coefficients, which
 * means \f$ Y_{00} \f$ must also be represented here as a SPHEREPACK
 * coefficient so that the control errors are expressed as SPHEREPACK
 * coefficients.
 *
 * \n
 * \note Since the \link domain::CoordinateMaps::TimeDependent::Shape Shape
 * \endlink map stores coefficients for all \f$ l \f$ but requires that the \f$
 * l=0 \f$ and \f$ l=1 \f$ coefficients be zero, the \f$ l=0 \f$ and \f$ l=1 \f$
 * modes of the control error are ecforced to be zero as well.
 */
template <::ah::ObjectLabel Horizon>
struct Shape : tt::ConformsTo<protocols::ControlError> {
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Computes the control error for shape control. This should not take any "
      "options."};

  void pup(PUP::er& /*p*/) {}

  template <typename Metavariables, typename... TupleTags>
  DataVector operator()(const Parallel::GlobalCache<Metavariables>& cache,
                        const double time,
                        const std::string& function_of_time_name,
                        const tuples::TaggedTuple<TupleTags...>& measurements) {
    const auto& domain = get<domain::Tags::Domain<3>>(cache);
    const auto& functions_of_time = get<domain::Tags::FunctionsOfTime>(cache);
    const DataVector lambda_lm_coefs =
        functions_of_time.at(function_of_time_name)->func(time)[0];
    const double lambda_00_coef =
        functions_of_time.at(detail::size_name<Horizon>())->func(time)[0][0];

    const auto& ah =
        get<control_system::QueueTags::Strahlkorper<Frame::Grid>>(measurements);
    const auto ah_coefs = ah.coefficients();

    ASSERT(lambda_lm_coefs.size() == ah_coefs.size(),
           "Number of coefficients for shape map '"
               << function_of_time_name << "' (" << lambda_lm_coefs.size()
               << ") does not math the number of coefficients in the AH "
                  "Strahlkorper ("
               << ah_coefs.size() << ").");

    const auto& excision_spheres = domain.excision_spheres();
    const double radius_excision_sphere_grid_frame =
        excision_spheres.at(detail::excision_sphere_name<Horizon>()).radius();

    // Y00 coefficient expressed as SPHEREPACK coefficient
    const double Y00_coef = detail::y00_coef();
    SpherepackIterator iter{ah.l_max(), ah.m_max()};
    const double relative_size_factor =
        (radius_excision_sphere_grid_frame / Y00_coef - lambda_00_coef) /
        ah_coefs[iter.set(0, 0)()];

    // The map parameters are in terms of SPHEREPACK coefficients (just like
    // strahlkorper coefficients), *not* spherical harmonic coefficients, thus
    // the control error for each l,m is in terms of SPHEREPACK coeffiecients
    // and no extra factors of sqrt(2/pi) are needed
    DataVector Q = -relative_size_factor * ah_coefs - lambda_lm_coefs;

    // Shape control is only for l > 1 so enfore that Q=0 for l=0,l=1
    Q[iter.set(0, 0)()] = 0.0;
    Q[iter.set(1, -1)()] = 0.0;
    Q[iter.set(1, 0)()] = 0.0;
    Q[iter.set(1, 1)()] = 0.0;

    return Q;
  }
};
}  // namespace ControlErrors
}  // namespace control_system
