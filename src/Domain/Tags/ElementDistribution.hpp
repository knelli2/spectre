// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>

#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementDistribution.hpp"
#include "Options/Auto.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

namespace domain {
namespace {
struct RoundRobin {};
}  // namespace

namespace OptionTags {
/// \ingroup OptionTagsGroup
/// \ingroup ComputationalDomainGroup
struct ElementDistribution {
  using type = Options::Auto<ElementWeight, RoundRobin>;
  static constexpr Options::String help = {
      "Weighting pattern to use for ZCurve element distribution. Specify "
      "RoundRobin to just place each element on the next core."};
};
}  // namespace OptionTags

namespace Tags {
struct ElementDistribution : db::SimpleTag {
  using type = std::optional<ElementWeight>;
  using option_tags = tmpl::list<OptionTags::ElementDistribution>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& element_distribution) {
    return element_distribution;
  }
};
}  // namespace Tags
}  // namespace domain
