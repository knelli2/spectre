// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/AsAccess.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Metafunctions.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp2 {
/*!
 * \brief Class that wraps a `db::Access` and allows it to be pupped.
 *
 * \details If we made the `db::Access` class pupable, then we'd have to
 * register every single `db::DataBox` class in our executable, which would be a
 * lot of extra work. This abstract base class and its derived classes get
 * around this issue by storing the typed `db::DataBox` in the derived class and
 * using the `AccessWrapper::pup_typed_box` virtual function to actually do the
 * pupping.
 */
struct AccessWrapper : public PUP::able {
  AccessWrapper() = default;
  WRAPPED_PUPable_abstract(AccessWrapper);

  virtual db::Access& access() = 0;

  void pup(PUP::er& p) { this->pup_typed_box(p); }

 private:
  virtual void pup_typed_box(PUP::er& p) = 0;
};

template <typename Target, size_t Dim>
struct TargetAccessWrapper : public AccessWrapper {
 private:
  using BoxType = intrp2::metafunctions::create_box_type<Target, Dim>;

 public:
  /// \cond
  TargetAccessWrapper() = default;
  explicit TargetAccessWrapper(CkMigrateMessage* /*unused*/) {}
  WRAPPED_PUPable_decl_template(TargetAccessWrapper);  // NOLINT
  /// \endcond

  TargetAccessWrapper(BoxType box) : box_(std::move(box)) {}

  db::Access& access() override { return db::as_access(box_); }

  void pup(PUP::er& p) override { p | box_; }

 private:
  void pup_typed_box(PUP::er& p) override { p | box_; }

  BoxType box_{};
};

/// \cond
template <typename Target, size_t Dim>
PUP::able::PUP_ID TargetAccessWrapper<Target, Dim>::my_PUP_ID = 0;  // NOLINT
/// \endcond
}  // namespace intrp2
