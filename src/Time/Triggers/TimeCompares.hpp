// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <pup.h>

#include "Options/Comparator.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct Time;
}  // namespace Tags
/// \endcond

namespace Triggers {
/// \ingroup EventsAndTriggersGroup
/// \ingroup TimeGroup
/// Trigger based on a comparison with the time.
class TimeCompares : public Trigger {
 public:
  /// \cond
  TimeCompares() = default;
  explicit TimeCompares(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(TimeCompares);  // NOLINT
  /// \endcond

  struct Comparison {
    using type = Options::Comparator;
    inline const static std::string help {"Comparison type"};
  };

  struct Value {
    using type = double;
    inline const static std::string help {"Value to compare to"};
  };

  using options = tmpl::list<Comparison, Value>;
  inline const static std::string help{
      "Trigger based on a comparison with the time."};

  explicit TimeCompares(const Options::Comparator comparator, const double time)
      : comparator_(comparator), time_(time) {}

  using argument_tags = tmpl::list<Tags::Time>;

  bool operator()(const double& time) const { return comparator_(time, time_); }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | comparator_;
    p | time_;
  }

 private:
  Options::Comparator comparator_{};
  double time_{};
};
}  // namespace Triggers
