// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <iomanip>
#include <sstream>
#include <string>

#include "ControlSystem/Tags/SystemTags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "Evolution/DiscontinuousGalerkin/InboxTags.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/OutputInbox.hpp"
#include "Parallel/Printf.hpp"
#include "Time/Actions/ChangeSlabSize.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace gh::deadlock {
struct PrintControlSystemCurrentMeasurement {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/) {
    const int current_measurement =
        db::get<control_system::Tags::CurrentNumberOfMeasurements>(box);

    Parallel::printf("%s: Current measurement = %d\n",
                     pretty_type::name<ParallelComponent>(),
                     current_measurement);
  }
};

struct PrintElementDeadlockInfo {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index) {
    auto& local_object =
        *Parallel::local(Parallel::get_parallel_component<ParallelComponent>(
            cache)[array_index]);

    const bool terminated = local_object.get_terminate();

    std::stringstream ss{};
    ss << std::scientific << std::setprecision(16);
    ss << "Element " << array_index
       << (terminated ? "terminated" : "did NOT terminate") << " at time "
       << db::get<::Tags::Time>(box) << ".";

    // Only print stuff if this element didn't terminate properly
    if (not terminated) {
      const std::string& next_action =
          local_object.deadlock_analysis_next_iterable_action();
      ss << " Next action: " << next_action << "\n";

      // Start with inboxes
      ss << " Inboxes:\n";

      const auto& inboxes = local_object.get_inboxes();

      const std::string mortar_inbox = Parallel::output_inbox<
          evolution::dg::Tags::BoundaryCorrectionAndGhostCellsInbox<3>>(inboxes,
                                                                        2_st);
      const std::string slab_size_inbox =
          Parallel::output_inbox<ChangeSlabSize_detail::NewSlabSizeInbox>(
              inboxes, 2_st);
      ss << mortar_inbox;
      ss << slab_size_inbox;

      ss << " Mortars:\n";

      const auto& mortar_next_temporal_id =
          db::get<evolution::dg::Tags::MortarNextTemporalId<3>>(box);

      ss << "  MortarNextTemporalId\n";
      for (const auto& [key, next_id] : mortar_next_temporal_id) {
        ss << "    Key: " << key << ", next time: " << next_id.substep_time()
           << "\n";
      }

      if constexpr (Metavariables::local_time_stepping) {
        const auto& mortar_data_history =
            db::get<evolution::dg::Tags::MortarDataHistory<
                3, typename db::add_tag_prefix<
                       ::Tags::dt,
                       typename Metavariables::system::variables_tag>::type>>(
                box);
        ss << "  MortarDataHistory:\n";

        for (const auto& [key, history] : mortar_data_history) {
          ss << "   Key: " << key << ", history:\n";
          ss << history.pretty_print_without_data(4_st);
        }
      } else {
        const auto& mortar_data =
            db::get<evolution::dg::Tags::MortarData<3>>(box);
        ss << "  MortarData:\n";

        for (const auto& [key, single_mortar_data] : mortar_data) {
          ss << "   Key: " << key << ", mortar data:\n";
          ss << single_mortar_data.pretty_print_current_buffer_no_data(4_st);
        }
      }
    } else {
      ss << "\n";
    }

    Parallel::printf(ss.str());
  }
};
}  // namespace gh::deadlock
