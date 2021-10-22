# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np

initial_separation = 15.366
initial_velocity = np.array([0.1, -0.2, 0.3])


def separation(time, init_sep):
    return (init_sep**4 - 12.8 * time)**0.25


def orbital_frequency(time, init_sep):
    return separation(time, init_sep)**-1.5


def positions1(time):
    init_sep = initial_separation
    return [
        0.5 * separation(time, init_sep) *
        np.cos(orbital_frequency(time, init_sep) * time) +
        initial_velocity[0] * time, 0.5 * separation(time, init_sep) *
        np.sin(orbital_frequency(time, init_sep) * time) +
        initial_velocity[1] * time, initial_velocity[2] * time
    ]


def positions2(time):
    init_sep = initial_separation
    return [
        -0.5 * separation(time, init_sep) *
        np.cos(orbital_frequency(time, init_sep) * time) +
        initial_velocity[0] * time, -0.5 * separation(time, init_sep) *
        np.sin(orbital_frequency(time, init_sep) * time) +
        initial_velocity[1] * time, initial_velocity[2] * time
    ]
