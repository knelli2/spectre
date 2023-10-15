// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <pup.h>
#include <pup_stl.h>
#include <string>
#include <utility>

#include "Options/Context.hpp"
#include "Options/Options.hpp"
#include "Options/ParseError.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Trigger.hpp"
#include "Time/Time.hpp"
#include "Time/TimeSequence.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Tags {
struct Time;
struct TimeStep;
}  // namespace Tags
/// \endcond

namespace Triggers {
namespace NearTimes_enums {
enum class Unit { Time, Slab, Step };
enum class Direction { Before, After, Both };
}  // namespace NearTimes_enums

/// \ingroup EventsAndTriggersGroup
/// \ingroup TimeGroup
/// Trigger in intervals surrounding particular times.
///
/// When using adaptive time stepping, intervals specified in terms of
/// slabs or steps are approximate.
///
/// \see Times
class NearTimes : public Trigger {
 public:
  /// \cond
  NearTimes() = default;
  explicit NearTimes(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(NearTimes);  // NOLINT
  /// \endcond

  using Unit = NearTimes_enums::Unit;
  using Direction = NearTimes_enums::Direction;

  struct OptionTags {
    struct Times {
      using type = std::unique_ptr<TimeSequence<double>>;
      inline const static std::string help {"Times to trigger at"};
    };

    struct Range {
      using type = double;
      static type lower_bound() { return 0.0; }
      inline const static std::string help
          {"Maximum time difference to trigger at"};
    };

    struct Unit {
      using type = NearTimes::Unit;
      inline const static std::string help
          {"Interpret Range as 'Time', 'Step's, or 'Slab's"};
    };

    struct Direction {
      using type = NearTimes::Direction;
      inline const static std::string help
          {"Trigger 'Before', 'After', or 'Both' from the times"};
    };
  };

  inline const static std::string help
      {"Trigger in intervals surrounding particular times."};
  using options =
      tmpl::list<typename OptionTags::Times, typename OptionTags::Range,
                 typename OptionTags::Unit, typename OptionTags::Direction>;

  NearTimes(std::unique_ptr<TimeSequence<double>> times, const double range,
            const Unit unit, const Direction direction)
      : times_(std::move(times)),
        range_(range),
        unit_(unit),
        direction_(direction) {}

  using argument_tags = tmpl::list<Tags::Time, Tags::TimeStep>;

  bool operator()(const double now, const TimeDelta& time_step) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  std::unique_ptr<TimeSequence<double>> times_{};
  double range_ = std::numeric_limits<double>::signaling_NaN();
  Unit unit_{};
  Direction direction_{};
};
}  // namespace Triggers

template <>
struct Options::create_from_yaml<Triggers::NearTimes_enums::Unit> {
  using type = Triggers::NearTimes_enums::Unit;
  template <typename Metavariables>
  static type create(const Options::Option& options) {
    const auto unit = options.parse_as<std::string>();
    if (unit == "Time") {
      return type::Time;
    } else if (unit == "Step") {
      return type::Step;
    } else if (unit == "Slab") {
      return type::Slab;
    } else {
      PARSE_ERROR(options.context(), "Unit must be 'Time', 'Step', or 'Slab'");
    }
  }
};

template <>
struct Options::create_from_yaml<
    typename Triggers::NearTimes_enums::Direction> {
  using type = Triggers::NearTimes_enums::Direction;
  template <typename Metavariables>
  static type create(const Options::Option& options) {
    const auto unit = options.parse_as<std::string>();
    if (unit == "Before") {
      return type::Before;
    } else if (unit == "After") {
      return type::After;
    } else if (unit == "Both") {
      return type::Both;
    } else {
      PARSE_ERROR(options.context(),
                  "Direction must be 'Before', 'After', or 'Both'");
    }
  }
};
