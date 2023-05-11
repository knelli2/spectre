
// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Cce {
template <bool RemoveFromInbox, typename Inbox, typename DbTags,
          typename Metavariables>
bool is_outer_boundary_element(
    const gsl::not_null<Inbox*> inbox, const db::DataBox<DbTags>& box,
    const Parallel::GlobalCache<Metavariables>& cache) {
  (void)inbox;  // Avoid unused variable warning
  const auto& element = db::get<domain::Tags::Element<3_st>>(box);
  const auto& external_boundaries = element.external_boundaries();
  // If we aren't on an external boundary clear the inbox and continue on.
  if (external_boundaries.size() == 0_st) {
    if constexpr (RemoveFromInbox) {
      inbox->clear();
    }
    return false;
  }

  const auto& domain = Parallel::get<domain::Tags::Domain<3>>(cache);
  const auto& excision_spheres = domain.excision_spheres();

  // If we are on the excision boundary, clear the inbox and continue on. We
  // are only concerned with the outer boundary.
  for (const auto& [name, excision_sphere] : excision_spheres) {
    (void)name;
    const auto& abutting_directions = excision_sphere.abutting_directions();
    for (const auto& element_direction : external_boundaries) {
      if (alg::any_of(abutting_directions,
                      [&element_direction](const auto& abutting_direction) {
                        return element_direction == abutting_direction.second;
                      })) {
        if constexpr (RemoveFromInbox) {
          inbox->clear();
        }
        return false;
      }
    }
  }

  return true;
}
}  // namespace Cce
