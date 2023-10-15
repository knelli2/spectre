// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <optional>
#include <pup.h>
#include <vector>

#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// Represents a sequence of times.
///
/// The template parameter \p T can either be `double` for a sequence
/// of simulation times or std::uint64_t for a sequence of slab
/// numbers.
template <typename T>
class TimeSequence : public PUP::able {
 protected:
  /// \cond HIDDEN_SYMBOLS
  TimeSequence() = default;
  TimeSequence(const TimeSequence&) = default;
  TimeSequence(TimeSequence&&) = default;
  TimeSequence& operator=(const TimeSequence&) = default;
  TimeSequence& operator=(TimeSequence&&) = default;
  /// \endcond

 public:
  ~TimeSequence() override = default;

  WRAPPED_PUPable_abstract(TimeSequence);  // NOLINT

  /// Returns the time in the sequence nearest to \p time, the time
  /// before that, and the time after, in numerical order.  These
  /// values allow a consumer to find the times in the sequence
  /// bracketing the given time and give enough additional information
  /// so that exact or roundoff matches can be handled however the
  /// consumer wishes.  Any time that does not exist (because the
  /// sequence terminates) is returned as an empty std::optional.  The
  /// central std::optional will be populated unless the sequence is
  /// empty.
  virtual std::array<std::optional<T>, 3> times_near(T time) const = 0;
};

/// \ingroup TimeGroup
///
/// Holds all the TimeSequences
namespace TimeSequences {
/// A sequence of evenly spaced times.
template <typename T>
class EvenlySpaced : public TimeSequence<T> {
 public:
  /// \cond
  EvenlySpaced() = default;
  explicit EvenlySpaced(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(EvenlySpaced);  // NOLINT
  /// \endcond

  struct Interval {
    inline const static std::string help {"Spacing between times"};
    using type = T;
    static constexpr T lower_bound() { return 0; }
  };

  struct Offset {
    inline const static std::string help {"Offset of sequence"};
    using type = T;
  };

  inline const static std::string help {"A sequence of evenly spaced times."};
  using options = tmpl::list<Interval, Offset>;

  explicit EvenlySpaced(T interval, T offset = 0,
                        const Options::Context& context = {});

  std::array<std::optional<T>, 3> times_near(T time) const override;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  // This is
  //   tmpl::conditional<std::is_integral_v<T>, std::make_signed<T>,
  //                     tmpl::identity<T>>::type
  // except avoiding instantiating std::make_signed<double> because,
  // depending on the standard version, that is either undefined
  // behavior or makes the program ill-formed.
  using SignedT = tmpl::apply<tmpl::apply<
      tmpl::if_<std::is_integral<tmpl::pin<T>>,
                tmpl::defer<tmpl::bind<std::make_signed_t, tmpl::pin<T>>>,
                tmpl::defer<tmpl::always<tmpl::pin<T>>>>>>;

  SignedT interval_{};
  SignedT offset_{};
};

/// An explicitly specified sequence of times.
template <typename T>
class Specified : public TimeSequence<T> {
 public:
  /// \cond
  Specified() = default;
  explicit Specified(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Specified);  // NOLINT
  /// \endcond

  struct Values {
    inline const static std::string help {"The times in the sequence"};
    using type = std::vector<T>;
  };

  inline const static std::string help
      {"An explicitly specified sequence of times."};
  using options = tmpl::list<Values>;

  explicit Specified(std::vector<T> values);

  std::array<std::optional<T>, 3> times_near(T time) const override;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  std::vector<T> values_;
};

template <typename T>
using all_time_sequences = tmpl::list<EvenlySpaced<T>, Specified<T>>;
}  // namespace TimeSequences
