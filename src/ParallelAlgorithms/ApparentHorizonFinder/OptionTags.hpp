// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <variant>
#include <vector>

#include "Domain/Structure/ObjectLabel.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/FastFlow.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Grid;
struct Distorted;
struct Inertial;
}  // namespace Frame
/// \endcond

namespace ah {
struct HorizonFinderOptions {
  struct InitialGuess {
    struct LMax {
      using type = size_t;
      static constexpr Options::String help = {
          "Strahlkorper is expanded in Ylms up to l=LMax"};
    };
    struct Radius {
      using type = double;
      static constexpr Options::String help = {
          "Radius of spherical Strahlkorper"};
    };
    struct Center {
      using type = std::array<double, 3>;
      static constexpr Options::String help = {
          "Center of spherical Strahlkorper"};
    };
    using options = tmpl::list<LMax, Radius, Center>;

    size_t l_max{};
    double radius{};
    std::array<double, 3> center{};

    // Option tag stuff
    static constexpr Options::String help = {"Initial guess"};
    using type = InitialGuess;
  };

  struct HorizonName {
    static constexpr Options::String help = {"Name of the horizon."};
    using type = domain::ObjectLabel;
    static std::string name() { return "Name"; }
  };
  struct HorizonFrame {
    static constexpr Options::String help = {
        "Frame of horizon find. Either Grid, Distorted, or Inertial."};
    using type = std::string;
    static std::string name() { return "Frame"; }
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

  using options =
      tmpl::list<InitialGuess, HorizonName, HorizonFrame, FastFlow, Verbosity>;
  static constexpr Options::String help = {
      "Provide an initial guess for the apparent horizon surface\n"
      "(Strahlkorper) and apparent-horizon-finding-algorithm (FastFlow)\n"
      "options."};

  HorizonFinderOptions(InitialGuess initial_guess_in, domain::ObjectLabel label,
                       const std::string& frame, ::FastFlow fast_flow_in,
                       ::Verbosity verbosity_in,
                       const Options::Context& context = {});

  HorizonFinderOptions() = default;
  HorizonFinderOptions(const HorizonFinderOptions& /*rhs*/) = default;
  HorizonFinderOptions& operator=(const HorizonFinderOptions& /*rhs*/) = delete;
  HorizonFinderOptions(HorizonFinderOptions&& /*rhs*/) = default;
  HorizonFinderOptions& operator=(HorizonFinderOptions&& /*rhs*/) = default;
  ~HorizonFinderOptions() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  std::variant<ylm::Strahlkorper<Frame::Grid>,
               ylm::Strahlkorper<Frame::Distorted>,
               ylm::Strahlkorper<Frame::Inertial>>
      initial_guess{};
  domain::ObjectLabel object_label{};
  ::FastFlow fast_flow{};
  ::Verbosity verbosity{::Verbosity::Quiet};
};

bool operator==(const HorizonFinderOptions& lhs,
                const HorizonFinderOptions& rhs);
bool operator!=(const HorizonFinderOptions& lhs,
                const HorizonFinderOptions& rhs);

namespace OptionTags {
struct HorizonFinders {
  using type = Options::Auto<std::vector<HorizonFinderOptions>,
                             Options::AutoLabel::None>;
  static constexpr Options::String help{"Options for a control system."};
  static std::string name() { return "ApparentHorizons"; }
};
}  // namespace OptionTags
}  // namespace ah
