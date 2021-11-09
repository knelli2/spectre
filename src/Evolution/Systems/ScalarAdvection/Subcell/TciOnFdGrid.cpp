// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarAdvection/Subcell/TciOnFdGrid.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/DgSubcell/PerssonTci.hpp"
#include "Evolution/Systems/ScalarAdvection/Subcell/TciOptions.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace ScalarAdvection::subcell {
template <size_t Dim>
bool TciOnFdGrid<Dim>::apply(const Scalar<DataVector>& dg_u,
                             const Mesh<Dim>& dg_mesh,
                             const TciOptions& tci_options,
                             const double persson_exponent) {
  const double max_abs_u = max(abs(get(dg_u)));
  return (max_abs_u > tci_options.u_cutoff) and
         ::evolution::dg::subcell::persson_tci(dg_u, dg_mesh, persson_exponent);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data) template struct TciOnFdGrid<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2))

#undef INSTANTIATION

#undef DIM
}  // namespace ScalarAdvection::subcell
