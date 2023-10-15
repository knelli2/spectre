// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <pup.h>

#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Options/String.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace domain::BoundaryConditions {
class MarkAsNone {
 public:
  MarkAsNone() = default;
  MarkAsNone(MarkAsNone&&) = default;
  MarkAsNone& operator=(MarkAsNone&&) = default;
  MarkAsNone(const MarkAsNone&) = default;
  MarkAsNone& operator=(const MarkAsNone&) = default;
  virtual ~MarkAsNone() = 0;
};

/*!
 * \brief None boundary conditions.
 *
 * This boundary condition doesn't actually do anything, and gets pretty much
 * completely ignored by everything but the domain creator internals. The domain
 * creator internals can use `None` as a way of specifying "boundary conditions"
 * without a system. It can also be used in cases like the BinaryCompactObject
 * domain where there may be no excision boundaries, and so the excision
 * boundary condition must be `None` in that case so the domain creator can be
 * sure the domain is set in a consistent state.
 *
 * To use with a specific system add:
 *
 * \code
 *  domain::BoundaryConditions::None<your::system::BoundaryConditionBase>
 * \endcode
 *
 * to the list of creatable classes.
 *
 * \warning if you want an outflow-type boundary condition, you must implement
 * one, not use `None.
 */
template <typename SystemBoundaryConditionBaseClass>
struct None final : public SystemBoundaryConditionBaseClass, public MarkAsNone {
 public:
  using options = tmpl::list<>;
  inline const static std::string help{
      "None boundary condition. Used only during domain creation to ensure a "
      "consistent state to the domain."};
  static std::string name() { return "None"; }

  None() = default;
  None(None&&) = default;
  None& operator=(None&&) = default;
  None(const None&) = default;
  None& operator=(const None&) = default;
  ~None() override = default;

  explicit None(CkMigrateMessage* msg);

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, None);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override;

  void pup(PUP::er& p) override;
};

template <typename SystemBoundaryConditionBaseClass>
None<SystemBoundaryConditionBaseClass>::None(CkMigrateMessage* const msg)
    : SystemBoundaryConditionBaseClass(msg) {}

template <typename SystemBoundaryConditionBaseClass>
std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
None<SystemBoundaryConditionBaseClass>::get_clone() const {
  return std::make_unique<None>(*this);
}

template <typename SystemBoundaryConditionBaseClass>
void None<SystemBoundaryConditionBaseClass>::pup(PUP::er& p) {
  BoundaryCondition::pup(p);
}

/// \cond
template <typename SystemBoundaryConditionBaseClass>
// NOLINTNEXTLINE
PUP::able::PUP_ID None<SystemBoundaryConditionBaseClass>::my_PUP_ID = 0;
/// \endcond

/// Check if a boundary condition inherits from `MarkAsNone`, which
/// constitutes as it being marked as a none boundary condition.
bool is_none(const std::unique_ptr<BoundaryCondition>& boundary_condition);
}  // namespace domain::BoundaryConditions
