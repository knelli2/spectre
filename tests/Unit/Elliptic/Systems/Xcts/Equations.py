# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np
from numpy import pi


def hamiltonian_sources(
    conformal_energy_density,
    extrinsic_curvature_trace,
    conformal_factor_minus_one,
    conformal_matter_scale=0,
):
    return (
        conformal_factor_minus_one + 1.0
    ) ** 5 * extrinsic_curvature_trace**2 / 12.0 - 2.0 * pi * (
        conformal_factor_minus_one + 1.0
    ) ** (
        5 - conformal_matter_scale
    ) * conformal_energy_density


def hamiltonian_sources_conf(*args, **kwargs):
    return hamiltonian_sources(*args, conformal_matter_scale=6, **kwargs)


def linearized_hamiltonian_sources(
    conformal_energy_density,
    extrinsic_curvature_trace,
    conformal_factor_minus_one,
    conformal_factor_correction,
    conformal_matter_scale=0,
):
    return (
        5.0
        * (conformal_factor_minus_one + 1.0) ** 4
        * extrinsic_curvature_trace**2
        / 12.0
        - 2.0
        * pi
        * (5.0 - conformal_matter_scale)
        * (conformal_factor_minus_one + 1.0) ** (4 - conformal_matter_scale)
        * conformal_energy_density
    ) * conformal_factor_correction


def linearized_hamiltonian_sources_conf(*args, **kwargs):
    return linearized_hamiltonian_sources(
        *args, conformal_matter_scale=6, **kwargs
    )


def distortion_hamiltonian_sources(
    longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
    conformal_factor_minus_one,
):
    return (
        -1.0
        / 32.0
        * (conformal_factor_minus_one + 1.0) ** 5
        * longitudinal_shift_minus_dt_conformal_metric_over_lapse_square
    )


def linearized_distortion_hamiltonian_sources(
    longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
    conformal_factor_minus_one,
    conformal_factor_correction,
):
    return (
        -5.0
        / 32.0
        * (conformal_factor_minus_one + 1.0) ** 4
        * conformal_factor_correction
        * longitudinal_shift_minus_dt_conformal_metric_over_lapse_square
    )


def curved_hamiltonian_or_lapse_sources(
    conformal_ricci_scalar, field, constant
):
    return (field + constant) * conformal_ricci_scalar / 8.0


def lapse_sources(
    conformal_energy_density,
    conformal_stress_trace,
    extrinsic_curvature_trace,
    dt_extrinsic_curvature_trace,
    shift_dot_deriv_extrinsic_curvature_trace,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_matter_scale=0,
):
    return (lapse_times_conformal_factor_minus_one + 1.0) * (
        conformal_factor_minus_one + 1.0
    ) ** 4 * (
        5.0 / 12.0 * extrinsic_curvature_trace**2
        + 2.0
        * pi
        * (conformal_energy_density + 2.0 * conformal_stress_trace)
        / (conformal_factor_minus_one + 1.0) ** conformal_matter_scale
    ) + (
        conformal_factor_minus_one + 1.0
    ) ** 5 * (
        shift_dot_deriv_extrinsic_curvature_trace - dt_extrinsic_curvature_trace
    )


def lapse_sources_conf(*args, **kwargs):
    return lapse_sources(*args, conformal_matter_scale=6, **kwargs)


def linearized_lapse_sources(
    conformal_energy_density,
    conformal_stress_trace,
    extrinsic_curvature_trace,
    dt_extrinsic_curvature_trace,
    shift_dot_deriv_extrinsic_curvature_trace,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_correction,
    lapse_times_conformal_factor_correction,
    conformal_matter_scale=0,
):
    return (
        (
            4.0
            * (lapse_times_conformal_factor_minus_one + 1.0)
            * (conformal_factor_minus_one + 1.0) ** 3
            * conformal_factor_correction
            + (conformal_factor_minus_one + 1.0) ** 4
            * lapse_times_conformal_factor_correction
        )
        * 5.0
        / 12.0
        * extrinsic_curvature_trace**2
        + (
            (4.0 - conformal_matter_scale)
            * (lapse_times_conformal_factor_minus_one + 1.0)
            * (conformal_factor_minus_one + 1.0) ** (3 - conformal_matter_scale)
            * conformal_factor_correction
            + (conformal_factor_minus_one + 1.0) ** (4 - conformal_matter_scale)
            * lapse_times_conformal_factor_correction
        )
        * 2.0
        * pi
        * (conformal_energy_density + 2.0 * conformal_stress_trace)
        + 5.0
        * (conformal_factor_minus_one + 1.0) ** 4
        * conformal_factor_correction
        * (
            shift_dot_deriv_extrinsic_curvature_trace
            - dt_extrinsic_curvature_trace
        )
    )


def linearized_lapse_sources_conf(*args, **kwargs):
    return linearized_lapse_sources(*args, conformal_matter_scale=6, **kwargs)


def distortion_hamiltonian_sources_with_lapse(
    longitudinal_shift_minus_dt_conformal_metric_square,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
):
    return (
        -1.0
        / 32.0
        * (conformal_factor_minus_one + 1.0) ** 7
        / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
        * longitudinal_shift_minus_dt_conformal_metric_square
    )


def linearized_distortion_hamiltonian_sources_with_lapse(
    longitudinal_shift_minus_dt_conformal_metric_square,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_correction,
    lapse_times_conformal_factor_correction,
):
    return (
        -1.0
        / 32.0
        * (
            7.0
            * (conformal_factor_minus_one + 1.0) ** 6
            * conformal_factor_correction
            / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
            - 2.0
            * (conformal_factor_minus_one + 1.0) ** 7
            / (lapse_times_conformal_factor_minus_one + 1.0) ** 3
            * lapse_times_conformal_factor_correction
        )
        * longitudinal_shift_minus_dt_conformal_metric_square
    )


def distortion_lapse_sources(
    longitudinal_shift_minus_dt_conformal_metric_square,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
):
    return (
        7.0
        / 32.0
        * (conformal_factor_minus_one + 1.0) ** 6
        / (lapse_times_conformal_factor_minus_one + 1.0)
        * longitudinal_shift_minus_dt_conformal_metric_square
    )


def linearized_distortion_lapse_sources(
    longitudinal_shift_minus_dt_conformal_metric_square,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_correction,
    lapse_times_conformal_factor_correction,
):
    return (
        7.0
        / 32.0
        * (
            6.0
            * (conformal_factor_minus_one + 1.0) ** 5
            / (lapse_times_conformal_factor_minus_one + 1.0)
            * conformal_factor_correction
            - (conformal_factor_minus_one + 1.0) ** 6
            / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
            * lapse_times_conformal_factor_correction
        )
        * longitudinal_shift_minus_dt_conformal_metric_square
    )


def momentum_sources(
    conformal_momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_metric,
    inv_conformal_metric,
    minus_div_dt_conformal_metric,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_flux,
    lapse_times_conformal_factor_flux,
    longitudinal_shift,
    conformal_matter_scale=0,
):
    return (
        np.einsum(
            "ij,jk,k",
            longitudinal_shift,
            conformal_metric,
            lapse_times_conformal_factor_flux
            / (lapse_times_conformal_factor_minus_one + 1.0)
            - 7.0 * conformal_factor_flux / (conformal_factor_minus_one + 1.0),
        )
        - minus_div_dt_conformal_metric
        + 4.0
        / 3.0
        * (lapse_times_conformal_factor_minus_one + 1.0)
        / (conformal_factor_minus_one + 1.0)
        * np.einsum(
            "ij,j", inv_conformal_metric, extrinsic_curvature_trace_gradient
        )
        + 16.0
        * pi
        * (lapse_times_conformal_factor_minus_one + 1.0)
        * (conformal_factor_minus_one + 1.0) ** (3 - conformal_matter_scale)
        * conformal_momentum_density
    )


def momentum_sources_conf(*args, **kwargs):
    return momentum_sources(*args, conformal_matter_scale=6, **kwargs)


def flat_cartesian_momentum_sources(
    conformal_momentum_density,
    extrinsic_curvature_trace_gradient,
    *args,
    **kwargs,
):
    return momentum_sources(
        conformal_momentum_density,
        extrinsic_curvature_trace_gradient,
        np.identity(3),
        np.identity(3),
        *args,
        **kwargs,
    )


def flat_cartesian_momentum_sources_conf(*args, **kwargs):
    return flat_cartesian_momentum_sources(
        *args, conformal_matter_scale=6, **kwargs
    )


def linearized_momentum_sources(
    conformal_momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_metric,
    inv_conformal_metric,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_flux,
    lapse_times_conformal_factor_flux,
    longitudinal_shift,
    conformal_factor_correction,
    lapse_times_conformal_factor_correction,
    shift_correction,
    conformal_factor_flux_correction,
    lapse_times_conformal_factor_flux_correction,
    longitudinal_shift_correction,
    conformal_matter_scale=0,
):
    return (
        np.einsum(
            "ij,jk,k",
            longitudinal_shift,
            conformal_metric,
            (
                lapse_times_conformal_factor_flux_correction
                / (lapse_times_conformal_factor_minus_one + 1.0)
                - lapse_times_conformal_factor_flux
                / ((lapse_times_conformal_factor_minus_one + 1.0) ** 2)
                * lapse_times_conformal_factor_correction
            )
            - 7.0
            * (
                conformal_factor_flux_correction
                / (conformal_factor_minus_one + 1.0)
                - conformal_factor_flux
                / (conformal_factor_minus_one + 1.0) ** 2
                * conformal_factor_correction
            ),
        )
        + np.einsum(
            "ij,jk,k",
            longitudinal_shift_correction,
            conformal_metric,
            lapse_times_conformal_factor_flux
            / (lapse_times_conformal_factor_minus_one + 1.0)
            - 7.0 * conformal_factor_flux / (conformal_factor_minus_one + 1.0),
        )
        + 4.0
        / 3.0
        * (
            lapse_times_conformal_factor_correction
            / (conformal_factor_minus_one + 1.0)
            - (lapse_times_conformal_factor_minus_one + 1.0)
            / (conformal_factor_minus_one + 1.0) ** 2
            * conformal_factor_correction
        )
        * np.einsum(
            "ij,j", inv_conformal_metric, extrinsic_curvature_trace_gradient
        )
        + 16.0
        * pi
        * (
            (3.0 - conformal_matter_scale)
            * (lapse_times_conformal_factor_minus_one + 1.0)
            * (conformal_factor_minus_one + 1.0) ** (2 - conformal_matter_scale)
            * conformal_factor_correction
            + (conformal_factor_minus_one + 1.0) ** (3 - conformal_matter_scale)
            * lapse_times_conformal_factor_correction
        )
        * conformal_momentum_density
    )


def linearized_momentum_sources_conf(*args, **kwargs):
    return linearized_momentum_sources(
        *args, conformal_matter_scale=6, **kwargs
    )


def flat_cartesian_linearized_momentum_sources(
    momentum_density, extrinsic_curvature_trace_gradient, *args, **kwargs
):
    return linearized_momentum_sources(
        momentum_density,
        extrinsic_curvature_trace_gradient,
        np.identity(3),
        np.identity(3),
        *args,
        **kwargs,
    )


def flat_cartesian_linearized_momentum_sources_conf(*args, **kwargs):
    return flat_cartesian_linearized_momentum_sources(
        *args, conformal_matter_scale=6, **kwargs
    )


def distortion_hamiltonian_sources_full(
    momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_metric,
    inv_conformal_metric,
    minus_div_dt_conformal_metric,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_flux,
    lapse_times_conformal_factor_flux,
    longitudinal_shift,
):
    return (
        -1.0
        / 32.0
        * (conformal_factor_minus_one + 1.0) ** 7
        / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
        * np.einsum(
            "ij,kl,ik,jl",
            conformal_metric,
            conformal_metric,
            longitudinal_shift,
            longitudinal_shift,
        )
    )


def flat_cartesian_distortion_hamiltonian_sources_full(
    momentum_density, extrinsic_curvature_trace_gradient, *args, **kwargs
):
    return distortion_hamiltonian_sources_full(
        momentum_density,
        extrinsic_curvature_trace_gradient,
        np.identity(3),
        np.identity(3),
        *args,
        **kwargs,
    )


def linearized_distortion_hamiltonian_sources_full(
    momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_metric,
    inv_conformal_metric,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_flux,
    lapse_times_conformal_factor_flux,
    longitudinal_shift,
    conformal_factor_correction,
    lapse_times_conformal_factor_correction,
    shift_correction,
    conformal_factor_flux_correction,
    lapse_times_conformal_factor_flux_correction,
    longitudinal_shift_correction,
):
    longitudinal_shift_square = np.einsum(
        "ij,kl,ik,jl",
        conformal_metric,
        conformal_metric,
        longitudinal_shift,
        longitudinal_shift,
    )
    longitudinal_shift_dot_correction = np.einsum(
        "ij,kl,ik,jl",
        conformal_metric,
        conformal_metric,
        longitudinal_shift,
        longitudinal_shift_correction,
    )
    return (
        -1.0
        / 32.0
        * (
            7.0
            * (conformal_factor_minus_one + 1.0) ** 6
            / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
            * longitudinal_shift_square
            * conformal_factor_correction
            - 2.0
            * (conformal_factor_minus_one + 1.0) ** 7
            / (lapse_times_conformal_factor_minus_one + 1.0) ** 3
            * longitudinal_shift_square
            * lapse_times_conformal_factor_correction
            + 2.0
            * (conformal_factor_minus_one + 1.0) ** 7
            / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
            * longitudinal_shift_dot_correction
        )
    )


def flat_cartesian_linearized_distortion_hamiltonian_sources_full(
    momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    *args,
    **kwargs,
):
    return linearized_distortion_hamiltonian_sources_full(
        momentum_density,
        extrinsic_curvature_trace_gradient,
        np.identity(3),
        np.identity(3),
        conformal_factor_minus_one,
        lapse_times_conformal_factor_minus_one,
        *args,
        **kwargs,
    )


def distortion_lapse_sources_with_shift(
    momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_metric,
    inv_conformal_metric,
    minus_div_dt_conformal_metric,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_flux,
    lapse_times_conformal_factor_flux,
    longitudinal_shift,
):
    return (
        7.0
        / 32.0
        * (conformal_factor_minus_one + 1.0) ** 6
        / (lapse_times_conformal_factor_minus_one + 1.0)
        * np.einsum(
            "ij,kl,ik,jl",
            conformal_metric,
            conformal_metric,
            longitudinal_shift,
            longitudinal_shift,
        )
    )


def flat_cartesian_distortion_lapse_sources_with_shift(
    momentum_density, extrinsic_curvature_trace_gradient, *args, **kwargs
):
    return distortion_lapse_sources_with_shift(
        momentum_density,
        extrinsic_curvature_trace_gradient,
        np.identity(3),
        np.identity(3),
        *args,
        **kwargs,
    )


def linearized_distortion_lapse_sources_with_shift(
    momentum_density,
    extrinsic_curvature_trace_gradient,
    conformal_metric,
    inv_conformal_metric,
    conformal_factor_minus_one,
    lapse_times_conformal_factor_minus_one,
    conformal_factor_flux,
    lapse_times_conformal_factor_flux,
    longitudinal_shift,
    conformal_factor_correction,
    lapse_times_conformal_factor_correction,
    shift_correction,
    conformal_factor_flux_correction,
    lapse_times_conformal_factor_flux_correction,
    longitudinal_shift_correction,
):
    longitudinal_shift_square = np.einsum(
        "ij,kl,ik,jl",
        conformal_metric,
        conformal_metric,
        longitudinal_shift,
        longitudinal_shift,
    )
    longitudinal_shift_dot_correction = np.einsum(
        "ij,kl,ik,jl",
        conformal_metric,
        conformal_metric,
        longitudinal_shift,
        longitudinal_shift_correction,
    )
    return 7.0 / 32.0 * (
        6.0
        * (conformal_factor_minus_one + 1.0) ** 5
        / (lapse_times_conformal_factor_minus_one + 1.0)
        * longitudinal_shift_square
        * conformal_factor_correction
        - (conformal_factor_minus_one + 1.0) ** 6
        / (lapse_times_conformal_factor_minus_one + 1.0) ** 2
        * longitudinal_shift_square
        * lapse_times_conformal_factor_correction
        + 2.0
        * (conformal_factor_minus_one + 1.0) ** 6
        / (lapse_times_conformal_factor_minus_one + 1.0)
        * longitudinal_shift_dot_correction
    ) + (conformal_factor_minus_one + 1.0) ** 5 * np.einsum(
        "i,i", shift_correction, extrinsic_curvature_trace_gradient
    )


def flat_cartesian_linearized_distortion_lapse_sources_with_shift(
    momentum_density, extrinsic_curvature_trace_gradient, *args, **kwargs
):
    return linearized_distortion_lapse_sources_with_shift(
        momentum_density,
        extrinsic_curvature_trace_gradient,
        np.identity(3),
        np.identity(3),
        *args,
        **kwargs,
    )
