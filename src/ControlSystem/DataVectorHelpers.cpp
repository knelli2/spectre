// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/DataVectorHelpers.hpp"

#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"

template <size_t N>
DataVector array_to_datavector(const std::array<double, N>& arr) {
  DataVector result{arr.size(), 0.0};
  for (size_t i = 0; i < N; i++) {
    result[i] = gsl::at(arr, i);
  }
  return result;
}

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

  result[0] = a[1] * b[2] - a[2] * b[1];
  result[1] = a[2] * b[0] - a[0] * b[2];
  result[2] = a[0] * b[1] - a[1] * b[0];

  return result;
}

template DataVector array_to_datavector(const std::array<double, 3>&);
