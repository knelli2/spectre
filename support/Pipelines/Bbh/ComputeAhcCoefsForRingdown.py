# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
from pathlib import Path
from typing import Optional, Union

import click
import numpy as np
import yaml
from rich.pretty import pretty_repr
from scipy.optimize import curve_fit as curve_fit

import spectre.IO.H5 as spectre_h5
from spectre.DataStructures import DataVector
from spectre.Domain import deserialize_functions_of_time
import spectre.Evolution.Ringdown as Ringdown

logger = logging.getLogger(__name__)


def compute_ahc_coefs_in_ringdown_distorted_frame(
    path_to_ah_h5,
    ahc_subfile_path,
    path_to_fot_h5,
    fot_subfile_path,
    path_to_output_h5,
    output_subfile_prefix,
    number_of_steps,
    which_obs_id,
    settling_timescale,
):
    output_subfile_ahc = output_subfile_prefix + "AhC_Ylm"
    output_subfile_dt_ahc = output_subfile_prefix + "dtAhC_Ylm"
    output_subfile_dt2_ahc = output_subfile_prefix + "dt2AhC_Ylm"

    ahc_times = []
    ahc_legend = ""
    ahc_center = []
    ahc_lmax = 0
    with spectre_h5.H5File(path_to_ah_h5, "r") as h5file:
        datfile = h5file.get_dat("/" + ahc_subfile_path)
        datfile_np = np.array(datfile.get_data())
        ahc_times = datfile_np[:, 0]
        ahc_legend = datfile.get_legend()
        ahc_center = [datfile_np[0][1], datfile_np[0][2], datfile_np[0][3]]
        ahc_lmax = int(datfile_np[0][4])

    # Transform AhC coefs to ringdown distorted frame
    with spectre_h5.H5File(path_to_fot_h5, "r") as h5file:
        volfile = h5file.get_vol(fot_subfile_path)
        obs_ids = volfile.list_observation_ids()
        fot_times = list(map(volfile.get_observation_value, obs_ids))
        functions_of_time = deserialize_functions_of_time(
            volfile.get_functions_of_time(obs_ids[which_obs_id])
        )

        # exp_func_and_2_derivs = [
        #     x[0] for x in functions_of_time['Expansion'].func_and_2_derivs(
        #         fot_times[which_obs_id])
        # ]
        exp_func_and_2_derivs = [1.0, 0.0, 0.0]

        exp_outer_bdry_func_and_2_derivs = [
            x[0]
            for x in functions_of_time[
                "ExpansionOuterBoundary"
            ].func_and_2_derivs(fot_times[which_obs_id])
        ]
        rot_func_and_2_derivs_tuple = functions_of_time[
            "Rotation"
        ].func_and_2_derivs(fot_times[which_obs_id])
        rot_func_and_2_derivs = [
            [coef for coef in x] for x in rot_func_and_2_derivs_tuple
        ]

        match_time = fot_times[which_obs_id]

        coefs_at_different_times = np.array(
            Ringdown.strahlkorper_coefs_in_ringdown_distorted_frame(
                path_to_ah_h5,
                ahc_subfile_path,
                ahc_times,
                number_of_steps,
                match_time,
                settling_timescale,
                exp_func_and_2_derivs,
                exp_outer_bdry_func_and_2_derivs,
                rot_func_and_2_derivs,
            )
        )

        # Print out coefficients for insertion into BBH domain
        print("Expansion: ", exp_func_and_2_derivs)
        print("ExpansionOutrBdry: ", exp_outer_bdry_func_and_2_derivs)
        print("Rotation: ", rot_func_and_2_derivs)
        print("Match time: ", match_time)
        print("Settling timescale: ", settling_timescale)
        print("Lmax: ", ahc_lmax)

    # Cubic fit transformed coefs to get first and second time derivatives
    def cubic(x, a, b, c, d):
        return a * x**3 + b * x**2 + c * x + d

    def dt_cubic(x, a, b, c, d):
        return 3 * a * x**2 + 2 * b * x + c

    def dt2_cubic(x, a, b, c, d):
        return 6 * a * x + 2 * b

    def fit_coefs(times, coefs):
        fits = []
        fit_ahc = []
        fit_dt_ahc = []
        fit_dt2_ahc = []
        for j in np.arange(0, coefs.shape[-1], 1):
            # if coef is not large enough, just set it to zero
            # Some coefficients are identically zero.
            # But I notice that my fits only recover the original coefficients
            # to ~1e-3 accuracy, so I decide to only fit coefficients bigger
            # than that, setting the rest to zero.
            if sum(np.abs(coefs[:, j])) < 1.0e-4:
                fits.append(np.zeros(4))
                fit_ahc.append(0.0)
                fit_dt_ahc.append(0.0)
                fit_dt2_ahc.append(0.0)
                continue
            fit = curve_fit(cubic, times, coefs[:, j])
            fits.append(fit[0])
            fit_ahc.append(cubic(match_time, *(fit[0])))
            fit_dt_ahc.append(dt_cubic(match_time, *(fit[0])))
            fit_dt2_ahc.append(dt2_cubic(match_time, *(fit[0])))
        return fit_ahc, fit_dt_ahc, fit_dt2_ahc

    # HACK: drop AhCs at times greater than the match time
    ahc_times_for_fit_list = []
    coefs_at_different_times_for_fit_list = []
    for i, time in enumerate(ahc_times[-number_of_steps:]):
        if time <= match_time:
            ahc_times_for_fit_list.append(time)
            coefs_at_different_times_for_fit_list.append(
                coefs_at_different_times[i]
            )
    ahc_times_for_fit = np.array(ahc_times_for_fit_list)
    coefs_at_different_times_for_fit = np.array(
        coefs_at_different_times_for_fit_list
    )
    print("AhC times available: " + str(ahc_times.shape))
    print("AhC times used: " + str(ahc_times_for_fit.shape))
    print("Coef times available: " + str(coefs_at_different_times.shape))
    print("Coef times used: " + str(coefs_at_different_times_for_fit.shape))

    fit_ahc_coefs, fit_ahc_dt_coefs, fit_ahc_dt2_coefs = fit_coefs(
        ahc_times[-number_of_steps:], coefs_at_different_times
    )

    # output coefs to H5
    # HACK: no translation, so inertial and distorted centers are the same,
    # i.e. are both the origin
    ahc_legend[1] = ahc_legend[1].replace("Inertial", "Distorted")
    ahc_legend[2] = ahc_legend[2].replace("Inertial", "Distorted")
    ahc_legend[3] = ahc_legend[3].replace("Inertial", "Distorted")

    fit_ahc_coefs_dv = -DataVector(fit_ahc_coefs)
    fit_ahc_dt_coefs_dv = -DataVector(fit_ahc_dt_coefs)
    fit_ahc_dt2_coefs_dv = -DataVector(fit_ahc_dt2_coefs)
    fit_ahc_coefs_to_write = Ringdown.wrap_fill_ylm_data(
        fit_ahc_coefs_dv, match_time, ahc_center, ahc_lmax
    )
    fit_ahc_dt_coefs_to_write = Ringdown.wrap_fill_ylm_data(
        fit_ahc_dt_coefs_dv, match_time, ahc_center, ahc_lmax
    )
    fit_ahc_dt2_coefs_to_write = Ringdown.wrap_fill_ylm_data(
        fit_ahc_dt2_coefs_dv, match_time, ahc_center, ahc_lmax
    )
    with spectre_h5.H5File(file_name=path_to_output_h5, mode="a") as h5file:
        ahc_datfile = h5file.insert_dat(
            path="/" + output_subfile_ahc, legend=ahc_legend, version=0
        )
        ahc_datfile.append(fit_ahc_coefs_to_write)

    with spectre_h5.H5File(file_name=path_to_output_h5, mode="a") as h5file:
        ahc_dt_datfile = h5file.insert_dat(
            path="/" + output_subfile_dt_ahc, legend=ahc_legend, version=0
        )
        ahc_dt_datfile.append(fit_ahc_dt_coefs_to_write)

    with spectre_h5.H5File(file_name=path_to_output_h5, mode="a") as h5file:
        ahc_dt2_datfile = h5file.insert_dat(
            path="/" + output_subfile_dt2_ahc, legend=ahc_legend, version=0
        )
        ahc_dt2_datfile.append(fit_ahc_dt2_coefs_to_write)


def compute_ahc_coefs_for_ringdown(
    path_to_ah_h5: Union[str, Path],
    ahc_subfile_path: str,
    path_to_fot_h5: Union[str, Path],
    fot_subfile_path: str,
    path_to_output_h5: Union[str, Path],
    output_subfile_prefix: str,
    number_of_steps: int,
    which_obs_id: int,
    settling_timescale: float,
):
    """Compute and write to disk ahc coefs in ringdown distorted frame."""
    logger.warning(
        "This code is still experimental! Ask Geoffrey if you have questions"
    )

    compute_ahc_coefs_in_ringdown_distorted_frame(
        str(path_to_ah_h5),
        ahc_subfile_path,
        str(path_to_fot_h5),
        fot_subfile_path,
        str(path_to_output_h5),
        output_subfile_prefix,
        number_of_steps,
        which_obs_id,
        settling_timescale,
    )


@click.command(
    name="compute-ahc-coefs-for-ringdown",
    help=compute_ahc_coefs_for_ringdown.__doc__,
)
@click.option(
    "-A",
    "--path_to_ah_h5",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        path_type=Path,
    ),
    help="Path to reduction file containing AhC coefs",
)
@click.option(
    "-a",
    "--ahc_subfile_path",
    type=str,
    help="Subfile path containing AhC coefs",
)
@click.option(
    "-F",
    "--path_to_fot_h5",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        path_type=Path,
    ),
    help="Path to volume data file containing functions of time",
)
@click.option(
    "-f",
    "--fot_subfile_path",
    type=str,
    help="Subfile path to functions of time",
)
@click.option(
    "-O",
    "--path_to_output_h5",
    type=click.Path(
        exists=False,
        file_okay=True,
        dir_okay=False,
        readable=True,
        path_type=Path,
    ),
    help="Path to output file",
)
@click.option(
    "-o", "--output_subfile_prefix", type=str, help="Output subfile prefix"
)
@click.option(
    "-n",
    "--number_of_steps",
    type=int,
    help="Number of steps from end to look for AhC data",
)
@click.option(
    "-w",
    "--which_obs_id",
    type=int,
    help="Which observation id to use for functions of time, matching time",
)
@click.option(
    "-s",
    "--settling_timescale",
    type=float,
    help="Damping timescale for settle to const",
)
def compute_ahc_coefs_for_ringdown_command(**kwargs):
    _rich_traceback_guard = True  # Hide traceback until here
    compute_ahc_coefs_for_ringdown(**kwargs)


if __name__ == "__main__":
    compute_ahc_coefs_for_ringdown_command(help_option_names=["-h", "--help"])
