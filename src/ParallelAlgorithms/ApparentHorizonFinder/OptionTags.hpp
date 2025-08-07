// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace ah {
/// Options for the adaptive horizon finder
struct AdaptivityOptions {
  static constexpr Options::String help = {
      "Options for an adaptive horizon finder."};
  using type = Options::Auto<AdaptivityOptions, Options::AutoLabel::None>;
  static std::string name() { return "Adaptivity"; }

  struct LBounds {
    static constexpr Options::String help = {
        "The minimum and maximum L used for the adaptive horizon finder. The "
        "max L must match the function of time L."};
    using type = std::array<size_t, 2>;
  };

  using options = tmpl::list<LBounds>;

  AdaptivityOptions() = default;

  AdaptivityOptions(std::array<size_t, 2> l_bounds_in,
                    const Options::Context& context = {});

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  std::array<size_t, 2> l_bounds{};
};

bool operator==(const AdaptivityOptions& lhs, const AdaptivityOptions& rhs);
bool operator!=(const AdaptivityOptions& lhs, const AdaptivityOptions& rhs);

/// Options for finding an apparent horizon.
template <typename Fr>
struct HorizonOptions {
 private:
  struct All {};

 public:
  /// See Strahlkorper for suboptions.
  struct InitialGuess {
    static constexpr Options::String help = {"Initial guess"};
    using type = ylm::Strahlkorper<Fr>;
  };
  /// See ::FastFlow for suboptions.
  struct FastFlow {
    static constexpr Options::String help = {"FastFlow options"};
    using type = ::FastFlow;
  };
  struct Verbosity {
    static constexpr Options::String help = {"Verbosity"};
    using type = ::Verbosity;
  };
  struct MaxComputeCoordsRetries {
    static constexpr Options::String help = {
        "Number of times to retry computing the coordinates of the horizon for "
        "each iteration. For the zeroth iteration, increases the 00 component "
        "by 50%. For subsequent iterations, two previous surfaces are averaged "
        "and that new surface is used."};
    using type = size_t;
  };
  struct BlocksForHorizonFind {
    static constexpr Options::String help = {
        "Volume data will be sent to the horizon finder from these block group "
        "names. Set to 'All' to send volume data from the entire domain."};
    using type = Options::Auto<std::vector<std::string>, All>;
  };
  using options =
      tmpl::list<InitialGuess, FastFlow, Verbosity, MaxComputeCoordsRetries,
                 BlocksForHorizonFind, AdaptivityOptions>;
  static constexpr Options::String help = {
      "Provide an initial guess for the apparent horizon surface\n"
      "(Strahlkorper) and apparent-horizon-finding-algorithm (FastFlow)\n"
      "options."};

  HorizonOptions(
      ylm::Strahlkorper<Fr> initial_guess_in, ::FastFlow fast_flow_in,
      ::Verbosity verbosity_in, size_t max_compute_coords_retries_in,
      std::optional<std::vector<std::string>> blocks_for_horizon_find_in,
      std::optional<AdaptivityOptions> adaptivity_in,
      const Options::Context& context = {});

  HorizonOptions() = default;
  HorizonOptions(const HorizonOptions& /*rhs*/) = default;
  HorizonOptions& operator=(const HorizonOptions& /*rhs*/) = delete;
  HorizonOptions(HorizonOptions&& /*rhs*/) = default;
  HorizonOptions& operator=(HorizonOptions&& /*rhs*/) = default;
  ~HorizonOptions() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  ylm::Strahlkorper<Fr> initial_guess{};
  ::FastFlow fast_flow;
  ::Verbosity verbosity{::Verbosity::Quiet};
  size_t max_compute_coords_retries{};
  std::optional<std::vector<std::string>> blocks_for_horizon_find;
  std::optional<AdaptivityOptions> adaptivity;
};

template <typename Fr>
bool operator==(const HorizonOptions<Fr>& lhs, const HorizonOptions<Fr>& rhs);
template <typename Fr>
bool operator!=(const HorizonOptions<Fr>& lhs, const HorizonOptions<Fr>& rhs);

namespace OptionTags {
struct ApparentHorizonGroup {
  static constexpr Options::String help{"Options for apparent horizon finders"};
  static std::string name() { return "ApparentHorizons"; }
};

template <typename HorizonMetavars>
struct ApparentHorizonOptions {
  using type = HorizonOptions<typename HorizonMetavars::frame>;
  static constexpr Options::String help{
      "Options for interpolation onto apparent horizon."};
  static std::string name() { return pretty_type::name<HorizonMetavars>(); }
  using group = ApparentHorizonGroup;
};
}  // namespace OptionTags
}  // namespace ah
