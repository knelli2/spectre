// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/WaveEquation/RegularSphericalWave.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <pup.h>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"

namespace ScalarWave::Solutions {

RegularSphericalWave::RegularSphericalWave(
    std::unique_ptr<MathFunction<1, Frame::Inertial>> profile)
    : profile_(std::move(profile)) {}

std::unique_ptr<evolution::initial_data::InitialData>
RegularSphericalWave::get_clone() const {
  return std::make_unique<RegularSphericalWave>(*this);
}

RegularSphericalWave::RegularSphericalWave(CkMigrateMessage* msg)
    : InitialData(msg) {}

RegularSphericalWave::RegularSphericalWave(const RegularSphericalWave& other)
    : evolution::initial_data::InitialData(other),
      profile_(other.profile_->get_clone()) {}

RegularSphericalWave& RegularSphericalWave::operator=(
    const RegularSphericalWave& other) {
  profile_ = other.profile_->get_clone();
  return *this;
}

constexpr double out_factor = 1.0;
constexpr double in_factor = 0.0;

tuples::TaggedTuple<Tags::Psi, Tags::Pi, Tags::Phi<3>>
RegularSphericalWave::variables(
    const tnsr::I<DataVector, 3>& x, double t,
    const tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3>> /*meta*/) const {
  const DataVector r = get(magnitude(x));
  // See class documentation for choice of cutoff
  const double r_cutoff = cbrt(std::numeric_limits<double>::epsilon());
  Scalar<DataVector> psi{r.size()};
  Scalar<DataVector> dpsi_dt{r.size()};
  tnsr::i<DataVector, 3> dpsi_dx{r.size()};
  for (size_t i = 0; i < r.size(); i++) {
    // Testing for r=0 here assumes a scale of order unity
    if (equal_within_roundoff(r[i], 0., r_cutoff, 1.)) {
      get(psi)[i] = 2. * profile_->first_deriv(-t);
      get(dpsi_dt)[i] = -2. * profile_->second_deriv(-t);
      for (size_t d = 0; d < 3; d++) {
        dpsi_dx.get(d)[i] = 0.;
      }
    } else {
      const auto F_out = out_factor * profile_->operator()(r[i] - t);
      const auto F_in = in_factor * profile_->operator()(-r[i] - t);
      const auto dF_out = out_factor * profile_->first_deriv(r[i] - t);
      const auto dF_in = in_factor * profile_->first_deriv(-r[i] - t);
      get(psi)[i] = (F_out - F_in) / r[i];
      get(dpsi_dt)[i] = (-dF_out + dF_in) / r[i];
      const double dpsi_dx_isotropic =
          (dF_out + dF_in - get(psi)[i]) / square(r[i]);
      for (size_t d = 0; d < 3; d++) {
        dpsi_dx.get(d)[i] = dpsi_dx_isotropic * x.get(d)[i];
      }
    }
  }
  tuples::TaggedTuple<Tags::Psi, Tags::Pi, Tags::Phi<3>> variables{
      std::move(psi), std::move(dpsi_dt), std::move(dpsi_dx)};
  get<Tags::Pi>(variables).get() *= -1.0;
  return variables;
}

tuples::TaggedTuple<::Tags::dt<Tags::Psi>, ::Tags::dt<Tags::Pi>,
                    ::Tags::dt<Tags::Phi<3>>>
RegularSphericalWave::variables(
    const tnsr::I<DataVector, 3>& x, double t,
    const tmpl::list<::Tags::dt<Tags::Psi>, ::Tags::dt<Tags::Pi>,
                     ::Tags::dt<Tags::Phi<3>>> /*meta*/) const {
  const DataVector r = get(magnitude(x));
  // See class documentation for choice of cutoff
  const double r_cutoff = cbrt(std::numeric_limits<double>::epsilon());
  Scalar<DataVector> dpsi_dt{r.size()};
  Scalar<DataVector> d2psi_dt2{r.size()};
  tnsr::i<DataVector, 3> d2psi_dtdx{r.size()};
  for (size_t i = 0; i < r.size(); i++) {
    // Testing for r=0 here assumes a scale of order unity
    if (equal_within_roundoff(r[i], 0., r_cutoff, 1.)) {
      get(dpsi_dt)[i] = -2. * profile_->second_deriv(-t);
      get(d2psi_dt2)[i] = 2. * profile_->third_deriv(-t);
      for (size_t d = 0; d < 3; d++) {
        d2psi_dtdx.get(d)[i] = 0.;
      }
    } else {
      const auto dF_out = out_factor * profile_->first_deriv(r[i] - t);
      const auto dF_in = in_factor * profile_->first_deriv(-r[i] - t);
      const auto d2F_out = out_factor * profile_->second_deriv(r[i] - t);
      const auto d2F_in = in_factor * profile_->second_deriv(-r[i] - t);
      get(dpsi_dt)[i] = (-dF_out + dF_in) / r[i];
      get(d2psi_dt2)[i] = (d2F_out - d2F_in) / r[i];
      const double d2psi_dtdx_isotropic =
          -(d2F_out + d2F_in + get(dpsi_dt)[i]) / square(r[i]);
      for (size_t d = 0; d < 3; d++) {
        d2psi_dtdx.get(d)[i] = d2psi_dtdx_isotropic * x.get(d)[i];
      }
    }
  }
  tuples::TaggedTuple<::Tags::dt<Tags::Psi>, ::Tags::dt<Tags::Pi>,
                      ::Tags::dt<Tags::Phi<3>>>
      dt_variables{std::move(dpsi_dt), std::move(d2psi_dt2),
                   std::move(d2psi_dtdx)};
  get<::Tags::dt<Tags::Pi>>(dt_variables).get() *= -1.0;
  return dt_variables;
}

bool operator==(const RegularSphericalWave& lhs,
                const RegularSphericalWave& rhs) {
  return *(lhs.profile_) == *(rhs.profile_);
}
bool operator!=(const RegularSphericalWave& lhs,
                const RegularSphericalWave& rhs) {
  return not(lhs == rhs);
}
void RegularSphericalWave::pup(PUP::er& p) {
  InitialData::pup(p);
  p | profile_;
}

PUP::able::PUP_ID RegularSphericalWave::my_PUP_ID = 0;
}  // namespace ScalarWave::Solutions
