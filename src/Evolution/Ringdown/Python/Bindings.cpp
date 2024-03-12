// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <pybind11/stl.h> if you want stl containers

#include "Domain/Creators/SphereTimeDependentMaps.hpp"
#include "Evolution/Ringdown/StrahlkorperCoefsInRingdownDistortedFrame.hpp"
#include "Utilities/ErrorHandling/SegfaultHandler.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace py = pybind11;

namespace evolution::Ringdown::py_bindings {  // NOLINT
void bind_strahlkorper_coefs_in_ringdown_distorted_frame(py::module& m) {
  tmpl::for_each<domain::creators::sphere::TimeDependentMapOptions::maps_list>(
      [](auto map_v) {
        using map = tmpl::type_from<decltype(map_v)>;
        register_classes_with_charm<map>();
      });
  m.def("strahlkorper_coefs_in_ringdown_distorted_frame",
        &evolution::Ringdown::strahlkorper_coefs_in_ringdown_distorted_frame,
        py::arg("path_to_horizons_h5"), py::arg("surface_subfile_name"),
        py::arg("requested_number_of_times_from_end"),
        py::arg("match_time"), py::arg("settling_timescale"),
        py::arg("exp_outer_bdry_func_and_2_derivs"),
        py::arg("rot_func_and_2_derivs"));
}
}  // namespace evolution::Ringdown::py_bindings

PYBIND11_MODULE(_Pybindings, m) {  // NOLINT
  enable_segfault_handler();
  // So I can return data vectors
  py::module_::import("spectre.DataStructures");
  evolution::Ringdown::py_bindings::
      bind_strahlkorper_coefs_in_ringdown_distorted_frame(m);
}
