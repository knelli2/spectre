// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace Convergence {

/*!
 * \brief Criteria that determine an iterative algorithm has converged
 *
 * \details Most criteria are based on a residual magnitude
 * \f$r_k\f$ after completion of an iteration \f$k\f$ (see, for instance, the
 * \ref LinearSolverGroup documentation, `LinearSolver::Tags::Residual` and
 * `LinearSolver::Tags::Magnitude`).
 *
 * The following criteria are implemented, ordered from highest to lowest
 * priority:
 *
 * - AbsoluteResidual: Matches if the residual has reached this magnitude.
 * - RelativeResidual: Matches if the residual has decreased by this factor,
 * relative to the start of the first iteration.
 * - MaxIterations: Matches if the number of iterations exceeds this limit.
 */
struct Criteria {
  inline const static std::string help
      {"The algorithm terminates when any of these criteria is matched."};

  struct MaxIterations {
    using type = size_t;
    inline const static std::string help {
        "The number of iterations exceeds this limit."};
    static type lower_bound() { return 0; }
  };

  struct AbsoluteResidual {
    using type = double;
    inline const static std::string help {
        "The residual has reached this magnitude."};
    static type lower_bound() { return 0.; }
  };

  struct RelativeResidual {
    using type = double;
    inline const static std::string help {
        "The residual has decreased by this factor."};
    static type lower_bound() { return 0.; }
    static type upper_bound() { return 1.; }
  };

  using options = tmpl::list<MaxIterations, AbsoluteResidual, RelativeResidual>;

  Criteria() = default;
  Criteria(size_t max_iterations_in, double absolute_residual_in,
           double relative_residual_in);

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  size_t max_iterations{};
  double absolute_residual{};
  double relative_residual{};
};

bool operator==(const Criteria& lhs, const Criteria& rhs);
bool operator!=(const Criteria& lhs, const Criteria& rhs);

}  // namespace Convergence
