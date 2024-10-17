// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <iomanip>
#include <ios>
#include <limits>
#include <sstream>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Printf/Printf.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolatedVars.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace deadlock {
struct PrintInterpolator {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(const db::DataBox<DbTags>& box,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index) {
    using target_tags = typename Metavariables::interpolation_target_tags;

    tmpl::for_each<target_tags>([&](const auto tag_v) {
      using target_tag = tmpl::type_from<decltype(tag_v)>;

      // Only need to print the sequential targets (aka horizons)
      if constexpr (target_tag::compute_target_points::is_sequential::value) {
        std::stringstream ss{};
        ss << std::setprecision(std::numeric_limits<double>::digits10 + 4)
           << std::scientific;

        const auto& holders =
            db::get<intrp::Tags::InterpolatedVarsHolders<Metavariables>>(box);
        const auto& vars_infos =
            get<intrp::Vars::HolderTag<target_tag, Metavariables>>(holders)
                .infos;
        const auto& num_elements =
            db::get<intrp::Tags::NumberOfElements<Metavariables::volume_dim>>(
                box);

        for (const auto& [temporal_id, info] : vars_infos) {
          ss.str("");
          ss << "Interpolator core " << array_index << ", "
             << pretty_type::name<target_tag>() << ", temporal id "
             << temporal_id << ": ";

          ss << "Iteration " << info.iteration << ", expecting data from "
             << num_elements.size() << " elements, but only received "
             << info.interpolation_is_done_for_these_elements.size()
             << ". Missing these elements:  ";

          std::unordered_set<ElementId<Metavariables::volume_dim>> difference{};
          for (const auto& element : num_elements) {
            if (info.interpolation_is_done_for_these_elements.count(element) ==
                0) {
              difference.insert(element);
            }
          }

          ss << difference;

          Parallel::printf("%s\n", ss.str());
        }
      }
    });

    (void)box;
    (void)cache;
    (void)array_index;
  }
};
}  // namespace deadlock
