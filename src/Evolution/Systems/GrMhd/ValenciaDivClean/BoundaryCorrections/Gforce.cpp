// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GrMhd/ValenciaDivClean/BoundaryCorrections/Gforce.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/NormalDotFlux.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace grmhd::ValenciaDivClean::BoundaryCorrections {
Gforce::Gforce(CkMigrateMessage* /*unused*/) {}

std::unique_ptr<BoundaryCorrection> Gforce::get_clone() const {
  return std::make_unique<Gforce>(*this);
}

void Gforce::pup(PUP::er& p) { BoundaryCorrection::pup(p); }

double Gforce::dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> packaged_tilde_d,
    const gsl::not_null<Scalar<DataVector>*> packaged_tilde_ye,
    const gsl::not_null<Scalar<DataVector>*> packaged_tilde_tau,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        packaged_tilde_s,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        packaged_tilde_b,
    const gsl::not_null<Scalar<DataVector>*> packaged_tilde_phi,
    const gsl::not_null<Scalar<DataVector>*> packaged_normal_dot_flux_tilde_d,
    const gsl::not_null<Scalar<DataVector>*> packaged_normal_dot_flux_tilde_ye,
    const gsl::not_null<Scalar<DataVector>*> packaged_normal_dot_flux_tilde_tau,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        packaged_normal_dot_flux_tilde_s,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        packaged_normal_dot_flux_tilde_b,
    const gsl::not_null<Scalar<DataVector>*> packaged_normal_dot_flux_tilde_phi,
    const gsl::not_null<Scalar<DataVector>*> packaged_abs_char_speed,
    const gsl::not_null<Scalar<DataVector>*> packaged_lapse,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        packaged_shift,
    const gsl::not_null<Scalar<DataVector>*> packaged_sqrt_det_spatial_metric,
    const gsl::not_null<tnsr::ii<DataVector, 3, Frame::Inertial>*>
        packaged_spatial_metric,
    const gsl::not_null<tnsr::II<DataVector, 3, Frame::Inertial>*>
        packaged_inv_spatial_metric,
    const gsl::not_null<Scalar<DataVector>*> packaged_pressure,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        packaged_spatial_velocity,
    const gsl::not_null<Scalar<DataVector>*> packaged_lorentz_factor,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        packaged_magnetic_field,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        packaged_normal_covector,

    const Scalar<DataVector>& tilde_d, const Scalar<DataVector>& tilde_ye,
    const Scalar<DataVector>& tilde_tau,
    const tnsr::i<DataVector, 3, Frame::Inertial>& tilde_s,
    const tnsr::I<DataVector, 3, Frame::Inertial>& tilde_b,
    const Scalar<DataVector>& tilde_phi,

    const tnsr::I<DataVector, 3, Frame::Inertial>& flux_tilde_d,
    const tnsr::I<DataVector, 3, Frame::Inertial>& flux_tilde_ye,
    const tnsr::I<DataVector, 3, Frame::Inertial>& flux_tilde_tau,
    const tnsr::Ij<DataVector, 3, Frame::Inertial>& flux_tilde_s,
    const tnsr::IJ<DataVector, 3, Frame::Inertial>& flux_tilde_b,
    const tnsr::I<DataVector, 3, Frame::Inertial>& flux_tilde_phi,

    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const Scalar<DataVector>& sqrt_det_spatial_metric,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inv_spatial_metric,

    const Scalar<DataVector>& pressure,
    const tnsr::I<DataVector, 3, Frame::Inertial>& spatial_velocity,
    const Scalar<DataVector>& lorentz_factor,
    const tnsr::I<DataVector, 3, Frame::Inertial>& magnetic_field,

    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& /*normal_vector*/,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
    /*mesh_velocity*/,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) {
  {
    // Compute max abs char speed
    Scalar<DataVector>& shift_dot_normal = *packaged_tilde_d;
    dot_product(make_not_null(&shift_dot_normal), shift, normal_covector);
    if (normal_dot_mesh_velocity.has_value()) {
      get(*packaged_abs_char_speed) =
          max(abs(-get(lapse) - get(shift_dot_normal) -
                  get(*normal_dot_mesh_velocity)),
              abs(get(lapse) - get(shift_dot_normal) -
                  get(*normal_dot_mesh_velocity)));
    } else {
      get(*packaged_abs_char_speed) =
          max(abs(-get(lapse) - get(shift_dot_normal)),
              abs(get(lapse) - get(shift_dot_normal)));
    }
  }

  *packaged_tilde_d = tilde_d;
  *packaged_tilde_ye = tilde_ye;
  *packaged_tilde_tau = tilde_tau;
  *packaged_tilde_s = tilde_s;
  *packaged_tilde_b = tilde_b;
  *packaged_tilde_phi = tilde_phi;

  normal_dot_flux(packaged_normal_dot_flux_tilde_d, normal_covector,
                  flux_tilde_d);
  normal_dot_flux(packaged_normal_dot_flux_tilde_ye, normal_covector,
                  flux_tilde_ye);
  normal_dot_flux(packaged_normal_dot_flux_tilde_tau, normal_covector,
                  flux_tilde_tau);
  normal_dot_flux(packaged_normal_dot_flux_tilde_s, normal_covector,
                  flux_tilde_s);
  normal_dot_flux(packaged_normal_dot_flux_tilde_b, normal_covector,
                  flux_tilde_b);
  normal_dot_flux(packaged_normal_dot_flux_tilde_phi, normal_covector,
                  flux_tilde_phi);

  *packaged_lapse = lapse;
  *packaged_shift = shift;
  *packaged_sqrt_det_spatial_metric = sqrt_det_spatial_metric;
  *packaged_spatial_metric = spatial_metric;
  *packaged_inv_spatial_metric = inv_spatial_metric;
  *packaged_pressure = pressure;
  *packaged_spatial_velocity = spatial_velocity;
  *packaged_lorentz_factor = lorentz_factor;
  *packaged_magnetic_field = magnetic_field;
  *packaged_normal_covector = normal_covector;

  return max(get(*packaged_abs_char_speed));
}

void Gforce::dg_boundary_terms(
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_tilde_d,
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_tilde_ye,
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_tilde_tau,
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*>
        boundary_correction_tilde_s,
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        boundary_correction_tilde_b,
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_tilde_phi,
    const Scalar<DataVector>& tilde_d_int,
    const Scalar<DataVector>& tilde_ye_int,
    const Scalar<DataVector>& tilde_tau_int,
    const tnsr::i<DataVector, 3, Frame::Inertial>& tilde_s_int,
    const tnsr::I<DataVector, 3, Frame::Inertial>& tilde_b_int,
    const Scalar<DataVector>& tilde_phi_int,
    const Scalar<DataVector>& normal_dot_flux_tilde_d_int,
    const Scalar<DataVector>& normal_dot_flux_tilde_ye_int,
    const Scalar<DataVector>& normal_dot_flux_tilde_tau_int,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_dot_flux_tilde_s_int,
    const tnsr::I<DataVector, 3, Frame::Inertial>& normal_dot_flux_tilde_b_int,
    const Scalar<DataVector>& normal_dot_flux_tilde_phi_int,
    const Scalar<DataVector>& abs_char_speed_int,
    const Scalar<DataVector>& lapse_int,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift_int,
    const Scalar<DataVector>& sqrt_det_spatial_metric_int,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric_int,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inv_spatial_metric_int,
    const Scalar<DataVector>& pressure_int,
    const tnsr::I<DataVector, 3, Frame::Inertial>& spatial_velocity_int,
    const Scalar<DataVector>& lorentz_factor_int,
    const tnsr::I<DataVector, 3, Frame::Inertial>& magnetic_field_int,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector_int,
    const Scalar<DataVector>& tilde_d_ext,
    const Scalar<DataVector>& tilde_ye_ext,
    const Scalar<DataVector>& tilde_tau_ext,
    const tnsr::i<DataVector, 3, Frame::Inertial>& tilde_s_ext,
    const tnsr::I<DataVector, 3, Frame::Inertial>& tilde_b_ext,
    const Scalar<DataVector>& tilde_phi_ext,
    const Scalar<DataVector>& normal_dot_flux_tilde_d_ext,
    const Scalar<DataVector>& normal_dot_flux_tilde_ye_ext,
    const Scalar<DataVector>& normal_dot_flux_tilde_tau_ext,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_dot_flux_tilde_s_ext,
    const tnsr::I<DataVector, 3, Frame::Inertial>& normal_dot_flux_tilde_b_ext,
    const Scalar<DataVector>& normal_dot_flux_tilde_phi_ext,
    const Scalar<DataVector>& abs_char_speed_ext,
    const Scalar<DataVector>& lapse_ext,
    const tnsr::I<DataVector, 3, Frame::Inertial>& shift_ext,
    const Scalar<DataVector>& sqrt_det_spatial_metric_ext,
    const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric_ext,
    const tnsr::II<DataVector, 3, Frame::Inertial>& inv_spatial_metric_ext,
    const Scalar<DataVector>& pressure_ext,
    const tnsr::I<DataVector, 3, Frame::Inertial>& spatial_velocity_ext,
    const Scalar<DataVector>& lorentz_factor_ext,
    const tnsr::I<DataVector, 3, Frame::Inertial>& magnetic_field_ext,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector_ext,
    const dg::Formulation dg_formulation) {
  if (dg_formulation == dg::Formulation::WeakInertial) {
    get(*boundary_correction_tilde_d) =
        0.5 * (get(normal_dot_flux_tilde_d_int) -
               get(normal_dot_flux_tilde_d_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_d_ext) - get(tilde_d_int));
    get(*boundary_correction_tilde_ye) =
        0.5 * (get(normal_dot_flux_tilde_ye_int) -
               get(normal_dot_flux_tilde_ye_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_ye_ext) - get(tilde_ye_int));
    get(*boundary_correction_tilde_tau) =
        0.5 * (get(normal_dot_flux_tilde_tau_int) -
               get(normal_dot_flux_tilde_tau_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_tau_ext) - get(tilde_tau_int));
    get(*boundary_correction_tilde_phi) =
        0.5 * (get(normal_dot_flux_tilde_phi_int) -
               get(normal_dot_flux_tilde_phi_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_phi_ext) - get(tilde_phi_int));

    for (size_t i = 0; i < 3; ++i) {
      boundary_correction_tilde_s->get(i) =
          0.5 * (normal_dot_flux_tilde_s_int.get(i) -
                 normal_dot_flux_tilde_s_ext.get(i)) -
          0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
              (tilde_s_ext.get(i) - tilde_s_int.get(i));
      boundary_correction_tilde_b->get(i) =
          0.5 * (normal_dot_flux_tilde_b_int.get(i) -
                 normal_dot_flux_tilde_b_ext.get(i)) -
          0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
              (tilde_b_ext.get(i) - tilde_b_int.get(i));
    }
  } else {
    get(*boundary_correction_tilde_d) =
        -0.5 * (get(normal_dot_flux_tilde_d_int) +
                get(normal_dot_flux_tilde_d_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_d_ext) - get(tilde_d_int));
    get(*boundary_correction_tilde_ye) =
        -0.5 * (get(normal_dot_flux_tilde_ye_int) +
                get(normal_dot_flux_tilde_ye_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_ye_ext) - get(tilde_ye_int));
    get(*boundary_correction_tilde_tau) =
        -0.5 * (get(normal_dot_flux_tilde_tau_int) +
                get(normal_dot_flux_tilde_tau_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_tau_ext) - get(tilde_tau_int));
    get(*boundary_correction_tilde_phi) =
        -0.5 * (get(normal_dot_flux_tilde_phi_int) +
                get(normal_dot_flux_tilde_phi_ext)) -
        0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
            (get(tilde_phi_ext) - get(tilde_phi_int));

    for (size_t i = 0; i < 3; ++i) {
      boundary_correction_tilde_s->get(i) =
          -0.5 * (normal_dot_flux_tilde_s_int.get(i) +
                  normal_dot_flux_tilde_s_ext.get(i)) -
          0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
              (tilde_s_ext.get(i) - tilde_s_int.get(i));
      boundary_correction_tilde_b->get(i) =
          -0.5 * (normal_dot_flux_tilde_b_int.get(i) +
                  normal_dot_flux_tilde_b_ext.get(i)) -
          0.5 * max(get(abs_char_speed_int), get(abs_char_speed_ext)) *
              (tilde_b_ext.get(i) - tilde_b_int.get(i));
    }
  }
}

bool operator==(const Gforce& /*lhs*/, const Gforce& /*rhs*/) { return true; }

bool operator!=(const Gforce& lhs, const Gforce& rhs) {
  return not(lhs == rhs);
}

// NOLINTNEXTLINE
PUP::able::PUP_ID Gforce::my_PUP_ID = 0;
}  // namespace grmhd::ValenciaDivClean::BoundaryCorrections
