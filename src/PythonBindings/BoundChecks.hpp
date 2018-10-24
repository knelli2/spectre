// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

#include "Utilities/PrettyType.hpp"

namespace py_bindings {
/*!
 * \ingroup PythonBindingsGroup
 * \brief Check if a vector-like object access is in bounds. Throws
 * std::runtime_error if it is not.
 */
template <typename T>
void bounds_check(const T& t, const size_t i) {
  if (i >= t.size()) {
    throw std::runtime_error{"Out of bounds access (" + std::to_string(i) +
                             ") into " + pretty_type::short_name<T>() +
                             " of size " + std::to_string(t.size())};
  }
}
}  // namespace py_bindings
