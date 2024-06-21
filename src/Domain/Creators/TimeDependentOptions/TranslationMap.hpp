// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedVariant.hpp"
#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

namespace domain::creators::time_dependent_options {
/*!
 * \brief Class to be used as an option for initializing translation map
 * coefficients.
 */
template <size_t Dim>
struct TranslationMapOptions {
  using type = Options::Auto<TranslationMapOptions, Options::AutoLabel::None>;
  static std::string name() { return "TranslationMap"; }
  static constexpr Options::String help = {
      "Options for a time-dependent translation of the coordinates. Specify "
      "'None' to not use this map."};

  struct Explicit {
    using type = std::array<std::array<double, Dim>, 3>;
    static constexpr Options::String help = "halp";
  };

  struct FromVolFile {
    using type = FromVolumeFile<names::Translation>;
    static constexpr Options::String help = "halp";
  };

  using Variant = variants::TaggedVariant<Explicit, FromVolFile>;

  struct InitialValues {
    using type = Variant;
    static constexpr Options::String help = {
        "Initial values for the translation map, its velocity and "
        "acceleration."};
  };

  using options = tmpl::list<InitialValues>;

  TranslationMapOptions() = default;
  TranslationMapOptions(Variant values_from_options);
  TranslationMapOptions(typename Explicit::type explicit_opt)
      : TranslationMapOptions(
            Variant{std::in_place_type<Explicit>, std::move(explicit_opt)}) {}
  TranslationMapOptions(typename FromVolFile::type from_vol_opt)
      : TranslationMapOptions(Variant{std::in_place_type<FromVolFile>,
                                      std::move(from_vol_opt)}) {}

  std::array<DataVector, 3> initial_values{};
};

/*!
 * \brief Given a `TranslationMapOptions`, return the function of time
 * coefficients that can be used to initialize a FunctionOfTime.
 */
template <size_t Dim>
std::array<DataVector, 3> initial_translation_func(
    const TranslationMapOptions<Dim>& translation_options);
}  // namespace domain::creators::time_dependent_options
