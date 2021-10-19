// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>

#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Options/Options.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/PrettyType.hpp"

namespace domain {
namespace creators {
namespace time_dependence {
/*!
 * \brief OptionTags for specific TimeDependences
 *
 * These are used to create Composition TimeDependences. E.g.
 * CompositionName:
 *   TimeDepOptionTag1:
 *     options...
 *   TimeDepOptionTag2:
 *     optons...
 */
namespace OptionTags {
/*!
 * OptionTag for all TimeDependences to be used in some kind of composition
 */
template <typename TimeDep>
struct TimeDependenceCompositionTag {
  static constexpr size_t mesh_dim = TimeDep::mesh_dim;
  static std::string name() { return TimeDep::name(); }
  using type = std::unique_ptr<TimeDependence<mesh_dim>>;
  static constexpr Options::String help = {
      "One of the maps in the composition."};
};
}  // namespace OptionTags
}  // namespace time_dependence
}  // namespace creators
}  // namespace domain
