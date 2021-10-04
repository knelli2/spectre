// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataVector.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

double dot(const DataVector& a, const DataVector& b);

DataVector cross(const DataVector& a, const DataVector& b);
