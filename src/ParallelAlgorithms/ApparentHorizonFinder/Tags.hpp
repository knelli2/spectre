// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <string>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/OptionTags.hpp"

/// \cond
class FastFlow;
/// \endcond

namespace ah::Tags {
struct FastFlow : db::SimpleTag {
  using type = ::FastFlow;
};

/// Base tag for whether or not to write the centers of the horizons to disk.
/// Most likely to be used in the `ObserveCenters` post horizon find callback
///
/// Other things can control whether the horizon centers are output by defining
/// their own simple tag from this base tag.
struct ObserveCentersBase : db::BaseTag {};

/// Simple tag for whether to write the centers of the horizons to disk.
/// Currently this tag is not creatable by options
struct ObserveCenters : ObserveCentersBase, db::SimpleTag {
  using type = bool;
};

struct HorizonFinders : db::SimpleTag {
  using type = std::optional<std::vector<HorizonFinderOptions>>;

  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::HorizonFinders>;

  static type create_from_options(const type& option) { return option; }
};
}  // namespace ah::Tags
