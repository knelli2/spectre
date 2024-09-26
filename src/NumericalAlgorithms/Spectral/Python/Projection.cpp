// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/Spectral/Python/Projection.hpp"

#include <cstddef>
#include <functional>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

#include "DataStructures/ApplyMatrices.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/Matrix.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Projection.hpp"

namespace py = pybind11;

namespace Spectral::py_bindings {
namespace {
template <size_t Dim>
void bind_projection_impl(py::module& m) {  // NOLINT
  m.def("p_projection_matrices", &p_projection_matrices<Dim>,
        py::arg("source_mesh"), py::arg("target_mesh"));
  m.def(
      "apply_matrices",
      [](const std::array<Matrix, Dim>& matrices, const DataVector& u,
         const Index<Dim>& extents) {
        return apply_matrices(matrices, u, extents);
      },
      py::arg("matrices"), py::arg("u"), py::arg("extents"));
}
}  // namespace

void bind_projection(py::module& m) {
  bind_projection_impl<1>(m);
  bind_projection_impl<2>(m);
  bind_projection_impl<3>(m);
}

}  // namespace Spectral::py_bindings
