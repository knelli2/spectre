// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"

template <size_t N>
DataVector array_to_datavector(const std::array<double, N>& arr);

double dot(const DataVector& a, const DataVector& b);

DataVector cross(const DataVector& a, const DataVector& b);
