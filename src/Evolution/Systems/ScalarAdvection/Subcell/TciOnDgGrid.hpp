// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/ScalarAdvection/Subcell/TciOptions.hpp"
#include "Evolution/Systems/ScalarAdvection/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
template <size_t Dim>
class Mesh;
/// \endcond

namespace ScalarAdvection::subcell {
/*!
 * \brief The troubled-cell indicator run on the DG grid to check if the
 * solution is admissible.
 *
 * Applies the Persson TCI to \f$U\f$ if its absolute value is greater than
 * `tci_options.u_cutoff`.
 */
template <size_t Dim>
struct TciOnDgGrid {
 public:
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<ScalarAdvection::Tags::U,
                                   ::domain::Tags::Mesh<Dim>, Tags::TciOptions>;

  static bool apply(const Scalar<DataVector>& dg_u, const Mesh<Dim>& dg_mesh,
                    const TciOptions& tci_options,
                    const double persson_exponent);
};
}  // namespace ScalarAdvection::subcell
