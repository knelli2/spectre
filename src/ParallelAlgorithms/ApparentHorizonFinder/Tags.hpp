// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <deque>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/OptionTags.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/Storage.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class FastFlow;
namespace Frame {
struct NoFrame;
}  // namespace Frame
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

struct VolumeDataAndCallbacks : db::SimpleTag {
  using type = std::map<Storage::NumberAndId, Storage::VolumeDataAndCallbacks>;
};

struct InterpolatedVars : db::simpleTag {
  using type = std::map<Storage::NumberAndId, Storage::InterpolatedVars>;
};

/// `temporal_id`s that we have already interpolated onto.
struct CompletedTemporalIds : db::SimpleTag {
  using type = std::set<Storage::NumberAndId>;
};

struct SurfaceAndPoints : db::simpleTag {
  using type =
      std::unordered_map<Storage::NumberAndId, Storage::SurfaceAndPoints>;
};

struct PreviousHorizons : db::SimpleTag {
  using type = std::deque<std::pair<double, ylm::Strahlkorper<Frame::NoFrame>>>;
};

struct Verbosity : db::SimpleTag {
  using type = ::Verbosity;
};
}  // namespace ah::Tags
