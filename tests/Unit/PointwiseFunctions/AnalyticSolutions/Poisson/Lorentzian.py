# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np

constant = 1.5


def field(x):
    return 1.0 / np.sqrt(1.0 + np.dot(x, x)) + constant


def field_gradient(x):
    return -np.asarray(x) / np.sqrt(1.0 + np.dot(x, x)) ** 3


def field_flux(x):
    return field_gradient(x)


def source(x):
    return 3.0 / np.sqrt(1.0 + np.dot(x, x)) ** 5
