// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Options/Auto.hpp"
#include "Options/Options.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"

namespace Parallel {
/// \ingroup ParallelGroup
/// Holds resource info for a single singleton.
template <typename ParallelComponent>
struct SingletonInfo {
  struct Proc {
    using type = Options::Auto<int>;
    static constexpr Options::String help = {
        "Proc to put singleton on. This can be determined automatically if "
        "desired."};
  };

  struct Exclusive {
    using type = bool;
    static constexpr Options::String help = {
        "Reserve this proc for this singleton. No DG elements will be placed "
        "on this proc."};
  };

  using options = tmpl::list<Proc, Exclusive>;
  static constexpr Options::String help = {
      "Resource options for all singletons."};

  SingletonInfo(std::optional<int> input_proc, const bool input_exclusive,
                const Options::Context& context = {})
      : proc_(input_proc ? *input_proc : -1), exclusive_(input_exclusive) {
    // If there is no value, we don't need to error so use 0 as a comparator
    // in both cases
    if (input_proc.value_or(0) < 0) {
      PARSE_ERROR(
          context,
          "Proc must be a non-negative integer. Please choose another proc.");
    }
    if (input_proc.value_or(0) > sys::number_of_procs() - 1) {
      PARSE_ERROR(context, "Proc (" << *input_proc
                                    << ") is beyond the last proc ("
                                    << sys::number_of_procs() - 1 << ").");
    }
  }

  SingletonInfo() = default;
  SingletonInfo(const SingletonInfo& /*rhs*/) = default;
  SingletonInfo& operator=(const SingletonInfo& /*rhs*/) = default;
  SingletonInfo(SingletonInfo&& /*rhs*/) = default;
  SingletonInfo& operator=(SingletonInfo&& /*rhs*/) = default;
  ~SingletonInfo() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | proc_;
    p | exclusive_;
  };

  int proc() const { return proc_; }
  bool exclusive() const { return exclusive_; }

 private:
  // Above, -1 is a sentinal value to represent that the proc should be chosen
  // automatically
  int proc_;
  bool exclusive_;
};

namespace OptionTags {
/// \ingroup ParallelGroup
/// Options group for singleton resource allocation
struct ResourceInfo {
  static std::string name() { return "ResourceInfo"; }
  static constexpr Options::String help = {"Options for allocating resources."};
};

/// \ingroup ParallelGroup
/// Option tag for a single singleton.
template <typename ParallelComponent>
struct SingletonInfo {
  using type = Parallel::SingletonInfo<ParallelComponent>;
  static constexpr Options::String help = {"Resource info for a singleton."};
  static std::string name() { return Options::name<ParallelComponent>(); }
  using group = ResourceInfo;
};

/// \ingroup ParallelGroup
/// Option tag to avoid placing Singleton and Array chares on proc 0
struct AvoidProc0 {
  using type = bool;
  static constexpr Options::String help = {
      "Whether to avoid placing DG elements on proc 0."};
  using group = ResourceInfo;
};
}  // namespace OptionTags

namespace Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup ParallelGroup
/// Tag that tells whether to avoid placing Singleton and Array chares on proc 0
struct AvoidProc0 : db::SimpleTag {
  using type = bool;
  using option_tags = tmpl::list<Parallel::OptionTags::AvoidProc0>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ParallelGroup
/// Tag that holds resource info about a singleton
template <typename ParallelComponent>
struct SingletonInfo : db::SimpleTag {
  using type = Parallel::SingletonInfo<ParallelComponent>;
  static std::string name() { return Options::name<ParallelComponent>(); }
  using option_tags =
      tmpl::list<Parallel::OptionTags::SingletonInfo<ParallelComponent>,
                 Parallel::OptionTags::AvoidProc0>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& singleton_info,
                                  const bool avoid_proc_0) {
    // Error if we are supposed to avoid proc 0, but a singleton has requested
    // to be placed on proc 0
    if (avoid_proc_0 and not singleton_info.proc()) {
      ERROR(
          "AvoidProc0 is true, but a singleton has requested to be placed on "
          "proc 0.");
    }
    return singleton_info;
  }
};
}  // namespace Tags
}  // namespace Parallel
