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
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace deadlock {
struct PrintInterpolationTarget {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(const db::DataBox<DbTags>& box,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/) {
    std::stringstream ss{};
    ss << std::setprecision(std::numeric_limits<double>::digits10 + 4)
       << std::scientific;

    using TargetTag = typename ParallelComponent::interpolation_target_tag;
    using TemporalId = typename TargetTag::temporal_id::type;

    const auto stream_points = [&](const TemporalId& temporal_id) {
      const auto& filled_indices =
          db::get<intrp::Tags::IndicesOfFilledInterpPoints<TemporalId>>(box);
      const auto& invalid_indices =
          db::get<intrp::Tags::IndicesOfInvalidInterpPoints<TemporalId>>(box);
      const auto& interpolated_vars =
          db::get<intrp::Tags::InterpolatedVars<TargetTag, TemporalId>>(box);

      const size_t expected_size =
          interpolated_vars.at(temporal_id).number_of_grid_points();
      const size_t filled_size = filled_indices.count(temporal_id) > 0
                                     ? filled_indices.at(temporal_id).size()
                                     : 0_st;
      const size_t invalid_size = invalid_indices.count(temporal_id) > 0
                                      ? invalid_indices.at(temporal_id).size()
                                      : 0_st;

      ss << "Total points expected = " << expected_size
         << ", valid points received = " << filled_size
         << ", invalid points received " << invalid_size << ". ";
    };

    if constexpr (TargetTag::compute_target_points::is_sequential::value) {
      ss << pretty_type::name<TargetTag>() << ", ";

      const auto& current_temporal_id =
          db::get<intrp::Tags::CurrentTemporalId<TemporalId>>(box);
      const auto& pending_temporal_ids =
          db::get<intrp::Tags::PendingTemporalIds<TemporalId>>(box);

      if (current_temporal_id.has_value()) {
        ss << "current temporal id " << current_temporal_id.value() << ", ";

        stream_points(current_temporal_id.value());
      } else {
        ss << "no current temporal id. ";
      }

      ss << "Pending ids " << pending_temporal_ids;

      Parallel::printf("%s\n", ss.str());
    } else {
      const auto& temporal_ids =
          db::get<intrp::Tags::TemporalIds<TemporalId>>(box);

      if (temporal_ids.empty()) {
        ss << pretty_type::name<TargetTag>() << ", No temporal ids.";
        Parallel::printf("%s\n", ss.str());
        return;
      }

      for (const auto& temporal_id : temporal_ids) {
        ss.str("");
        ss << pretty_type::name<TargetTag>() << ", temporal id " << temporal_id
           << ", ";

        stream_points(temporal_id);

        Parallel::printf("%s\n", ss.str());
      }
    }
  }
};
}  // namespace deadlock
