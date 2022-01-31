// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Tags.hpp"
#include "Evolution/DgSubcell/Tags/Inactive.hpp"
#include "Evolution/Systems/ScalarAdvection/Subcell/TciOptions.hpp"
#include "Evolution/Systems/ScalarAdvection/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t Dim>
class Mesh;
template <typename TagsList>
class Variables;
/// \endcond

namespace ScalarAdvection::subcell {
/*!
 * \brief The troubled-cell indicator run on DG initial data to see if we need
 * to switch to subcell.
 *
 * - Apply the Persson TCI to the scalar field \f$U\f$ if its magnitude on the
 * DG grid is greater than `tci_options.u_cutoff`.
 * - Apply the two-mesh relaxed discrete maximum principle TCI.
 */
template <size_t Dim>
struct DgInitialDataTci {
 private:
  template <typename Tag>
  using Inactive = evolution::dg::subcell::Tags::Inactive<Tag>;

 public:
  using argument_tags = tmpl::list<domain::Tags::Mesh<Dim>, Tags::TciOptions>;

  static bool apply(
      const Variables<tmpl::list<ScalarAdvection::Tags::U>>& dg_vars,
      const Variables<tmpl::list<Inactive<ScalarAdvection::Tags::U>>>&
          subcell_vars,
      double rdmp_delta0, double rdmp_epsilon, double persson_exponent,
      const Mesh<Dim>& dg_mesh, const TciOptions& tci_options);
};
}  // namespace ScalarAdvection::subcell
