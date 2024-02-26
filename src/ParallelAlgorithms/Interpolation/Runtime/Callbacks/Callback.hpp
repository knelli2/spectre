// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <pup.h>

#include "DataStructures/DataBox/Access.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/CallWithDynamicType.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp2::callbacks {
// The Callback base class is templated on the Target because different targets
// will have different sets of callbacks, and they need to be distinguishable in
// option creation.
template <typename Target>
class Callback : PUP::able {
  /// \cond
  Callback() = default;
  Callback(const Callback&) = default;
  Callback(Callback&&) = default;
  Callback& operator=(const Callback&) = default;
  Callback& operator=(Callback&&) = default;
  /// \endcond

  ~Callback() override = default;
  explicit Callback(CkMigrateMessage* msg) : PUP::able(msg) {}
  WRAPPED_PUPable_abstract(Callback);  // NOLINT

  template <typename Metavariables, typename TemporalId>
  void invoke(const db::Access& access,
              const Parallel::GlobalCache<Metavariables>& cache,
              const TemporalId& temporal_id) {
    using factory_classes =
        typename std::decay_t<Metavariables>::factory_creation::factory_classes;
    return call_with_dynamic_type<void,
                                  tmpl::at<factory_classes, Callback<Target>>>(
        this, [&](auto* const derived_callback) {
          derived_callback->apply(access, cache, temporal_id);
        });
  }

  virtual const std::unordered_set<std::string>& observables() const = 0;
};
}  // namespace intrp2::callbacks
