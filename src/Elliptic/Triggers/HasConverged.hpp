// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <pup.h>

#include "NumericalAlgorithms/Convergence/HasConverged.hpp"
#include "NumericalAlgorithms/Convergence/Tags.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace elliptic::Triggers {
/// \ingroup EventsAndTriggersGroup
/// Trigger when the solver identified by the `Label` has converged.
template <typename Label>
class HasConverged : public Trigger {
 public:
  /// \cond
  HasConverged() = default;
  explicit HasConverged(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(HasConverged);  // NOLINT
  /// \endcond

  using options = tmpl::list<>;
  inline const static std::string help{
      "Trigger when the solver has converged."};

  using argument_tags = tmpl::list<Convergence::Tags::HasConverged<Label>>;

  bool operator()(const Convergence::HasConverged& has_converged) const {
    return static_cast<bool>(has_converged);
  }
};

/// \cond
template <typename Label>
PUP::able::PUP_ID HasConverged<Label>::my_PUP_ID = 0;  // NOLINT
/// \endcond
}  // namespace elliptic::Triggers
