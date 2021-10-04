// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataVector.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

double dot(const DataVector& a, const DataVector& b) {
  ASSERT(a.size() == b.size(),
         "DataVectors must be the same size to dot together.");
  ASSERT(a.size() != 0,
         "DataVectors can't be size 0 if you want to dot them together.");
  double result = 0.0;
  for (size_t i = 0; i < a.size(); i++) {
    result += a[i] * b[i];
  }
  return result;
}

DataVector cross(const DataVector& a, const DataVector& b) {
  ASSERT(a.size() == 3, "Can only cross DataVectors that are size 3.");
  ASSERT(a.size() == b.size(),
         "DataVectors must be the same size to cross together.");
  DataVector result{3, 0.0};

  result[0] = a[0] * b[1] - a[1] * b[0];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[0] * b[1];

  return result;
}
