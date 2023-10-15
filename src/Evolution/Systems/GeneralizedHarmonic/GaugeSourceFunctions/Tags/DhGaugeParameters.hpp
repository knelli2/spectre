// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>

#include "Evolution/Systems/GeneralizedHarmonic/GaugeSourceFunctions/DhGaugeParameters.hpp"
#include "Options/String.hpp"

/// \cond
namespace gh::OptionTags {
struct Group;
}  // namespace gh::OptionTags
/// \endcond

namespace gh::gauges {
namespace OptionTags {
template <bool UseRollon>
struct DhGaugeParameters {
  using type = gh::gauges::DhGaugeParameters<UseRollon>;
  inline const static std::string help{
      "Parameters for initializing damped harmonic gauge."};
  using group = gh::OptionTags::Group;
};
}  // namespace OptionTags

namespace Tags {
/// \brief Input option tags for the generalized harmonic evolution system
template <bool UseRollon>
struct DhGaugeParameters : db::SimpleTag {
  using ParametersType = gh::gauges::DhGaugeParameters<UseRollon>;
  using type = ParametersType;
  using option_tags =
      tmpl::list<gh::gauges::OptionTags::DhGaugeParameters<UseRollon>>;

  static constexpr bool pass_metavariables = false;
  static ParametersType create_from_options(const ParametersType& parameters) {
    return parameters;
  }
};
}  // namespace Tags
}  // namespace gh::gauges
