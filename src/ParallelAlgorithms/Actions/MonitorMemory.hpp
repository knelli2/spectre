// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <numeric>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/MemoryMonitor/MemoryMonitor.hpp"
#include "Parallel/MemoryMonitor/Tags.hpp"
#include "Parallel/Serialize.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"

namespace mem_monitor {
/*!
 * \brief Simple action meant to be run on the MemoryMonitor component.
 *
 * \details This action collects the sizes of all the local branches of a group
 * or nodegroup component, computes the total memory usage on a node for each,
 * then writes it to disk. For groups, the proc with the maximum memory usage is
 * also reported along with the size on the proc.
 *
 * The columns in the dat file for a nodegroup when running on 3 nodes will be
 *
 * - %Time
 * - Size on node 0 (MB)
 * - Size on node 1 (MB)
 * - Size on node 2 (MB)
 * - Average size per node (MB)
 *
 * The columns in the dat file for a group when running on 3 nodes will be
 *
 * - %Time
 * - Size on node 0 (MB)
 * - Size on node 1 (MB)
 * - Size on node 2 (MB)
 * - Proc of max size
 * - Size on proc of max size (MB)
 * - Average size per node (MB)
 *
 * The dat file will be placed in the `/MemoryMonitors/` group in the reduction
 * file. The name of the dat file is the `pretty_type::name` of the component.
 */
template <typename ContributingComponent>
struct ContributeMemoryData {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time,
                    const int node_or_proc, const double size_in_MB) {
    using tag = Tags::MemoryHolder<ContributingComponent>;
    if constexpr (db::tag_is_retrievable_v<tag, db::DataBox<DbTags>>) {
      db::mutate<tag>(
          make_not_null(&box),
          [&cache, &time, &node_or_proc, &size_in_MB](
              gsl::not_null<
                  std::unordered_map<double, std::unordered_map<int, double>>*>
                  memory_holder) {
            // If this is the first reduction at this time, create a new
            // unordered map
            if (memory_holder->count(time) == 0) {
              memory_holder->emplace(time, std::unordered_map<int, double>{
                                               {node_or_proc, size_in_MB}});
            } else {
              // otherwise continue to add to the existing time
              memory_holder->at(time)[node_or_proc] = size_in_MB;
            }

            // If we have received data for every node/proc at a given
            // time, get all the data, write it to disk, then remove the current
            // time from the stored times as it's no longer needed

            auto& mem_monitor_proxy =
                Parallel::get_parallel_component<MemoryMonitor<Metavariables>>(
                    cache);

            constexpr bool is_group =
                Parallel::is_group_v<ContributingComponent>;

            const int num_nodes =
                Parallel::number_of_nodes(*Parallel::local(mem_monitor_proxy));
            const int num_procs =
                Parallel::number_of_procs(*Parallel::local(mem_monitor_proxy));
            const size_t expected_number =
                static_cast<size_t>(is_group ? num_procs : num_nodes);
            ASSERT(memory_holder->at(time).size() <= expected_number,
                   "ContributeMemoryData received more data that it was "
                   "expecting. Was expecting "
                       << expected_number << " calls but instead got "
                       << memory_holder->at(time).size());
            if (memory_holder->at(time).size() == expected_number) {
              // First column is always time
              std::vector<double> data_to_append{time};
              std::vector<std::string> legend{{"Time"}};

              // Append a column for each node, and keep track of cumulative
              // total. If we have proc data (from groups) do an additional loop
              // over the procs to get the total on that node and get the proc
              // of the maximum memory usage
              double avg_size_per_node = 0.0;
              double max_usage_on_proc = -std::numeric_limits<double>::max();
              int proc_of_max = 0;
              for (int node = 0; node < num_nodes; node++) {
                double size_on_node = 0.0;
                if (not is_group) {
                  size_on_node = memory_holder->at(time).at(node);
                } else {
                  const int first_proc = Parallel::first_proc_on_node(
                      node, *Parallel::local(mem_monitor_proxy));
                  const int procs_on_node = Parallel::procs_on_node(
                      node, *Parallel::local(mem_monitor_proxy));
                  const int last_proc = first_proc + procs_on_node;
                  for (int proc = first_proc; proc < last_proc; proc++) {
                    size_on_node += memory_holder->at(time).at(proc);
                    if (memory_holder->at(time).at(proc) > max_usage_on_proc) {
                      max_usage_on_proc = memory_holder->at(time).at(proc);
                      proc_of_max = proc;
                    }
                  }
                }

                data_to_append.push_back(size_on_node);
                avg_size_per_node += size_on_node;
                legend.emplace_back("Size on node " + get_output(node) +
                                    " (MB)");
              }

              // If we have proc data, write the proc with the maximum usage to
              // disk along with how much memory it's using
              if (is_group) {
                data_to_append.push_back(static_cast<double>(proc_of_max));
                data_to_append.push_back(max_usage_on_proc);
                legend.emplace_back("Proc of max size");
                legend.emplace_back("Size on proc of max size (MB)");
              }

              avg_size_per_node /= static_cast<double>(num_nodes);

              // Last column is average over all nodes
              data_to_append.push_back(avg_size_per_node);
              legend.emplace_back("Average size per node (MB)");

              auto& observer_writer_proxy = Parallel::get_parallel_component<
                  observers::ObserverWriter<Metavariables>>(cache);

              Parallel::threaded_action<
                  observers::ThreadedActions::WriteReductionDataRow>(
                  // Node 0 is always the writer
                  observer_writer_proxy[0],
                  subfile_name<ContributingComponent>(), legend,
                  std::make_tuple(data_to_append));

              // Clean up finished time
              auto finished_time_iter = memory_holder->find(time);
              memory_holder->erase(finished_time_iter);
            }
          });
    } else {
      (void)box;
      (void)cache;
      (void)time;
      (void)node_or_proc;
      (void)size_in_MB;
      ERROR(
          "Wrong DataBox. Expected the DataBox for the MemoryMonitor which has "
          "the mem_monitor::Tags::MemoryHolder tag.");
    }
  }
};

/*!
 * \brief Simple action meant to be used as a callback for
 * Parallel::contribute_to_reduction that writes the size of an Array parallel
 * component to disk.
 *
 * \details The columns in the dat file when running on 3 nodes will be
 *
 * - %Time
 * - Size on node 0 (MB)
 * - Size on node 1 (MB)
 * - Size on node 2 (MB)
 * - Average size per node (MB)
 *
 * The dat file will be placed in the `/MemoryMonitors/` group in the reduction
 * file. The name of the dat file is the `pretty_type::name` of the component.
 */
template <typename ArrayComponent>
struct ProcessArray {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time,
                    const std::vector<double>& size_per_node) {
    auto& observer_writer_proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);

    std::vector<std::string> legend{{"Time"}};
    for (size_t i = 0; i < size_per_node.size(); i++) {
      const std::string label = "Size on node " + get_output(i) + " (MB)";
      legend.emplace_back(label);
    }
    legend.emplace_back("Average size per node (MB)");

    const double avg_size =
        std::accumulate(size_per_node.begin(), size_per_node.end(), 0.0) /
        static_cast<double>(size_per_node.size());

    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        // Node 0 is always the writer
        observer_writer_proxy[0], subfile_name<ArrayComponent>(), legend,
        std::make_tuple(time, size_per_node, avg_size));
  }
};

/*!
 * \brief Simple action meant to be run on every branch of a Nodegroup that
 * computes the size of the local branch and reports that size to the
 * MemoryMonitor using the ContributeMemory simple action.
 */
struct ProcessNodegroup {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time) {
    auto& singleton_proxy =
        Parallel::get_parallel_component<MemoryMonitor<Metavariables>>(cache);

    auto& nodegroup_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const int my_node =
        Parallel::my_node(*Parallel::local_branch(nodegroup_proxy));

    const double size_in_bytes = static_cast<double>(
        size_of_object_in_bytes(*Parallel::local_branch(nodegroup_proxy)));

    const double size_in_MB = size_in_bytes / 1.0e6;

    // Note that we don't call Parallel::contribute_to_reduction here like we
    // do for charm Arrays in the MonitorMemory Event. This is because Charm
    // requires that all calls to contribute_to_reduction happen in the exact
    // same order every time, otherwise it is undefined behavior. However,
    // this simple action (ProcessNodegroup) is called on each branch of the
    // nodegroup. Because this is a simple action, the order that the branches
    // run the simple actions is completely random based on communication
    // patterns in charm. Thus, calling contribute_to_reduction here would
    // result in undefined behavior.
    Parallel::simple_action<ContributeMemoryData<ParallelComponent>>(
        singleton_proxy, time, my_node, size_in_MB);
  }
};

/*!
 * \brief Simple action meant to be run on every branch of a Group that
 * computes the size of the local branch and reports that size to the
 * MemoryMonitor using the ContributeMemory simple action.
 */
struct ProcessGroup {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time) {
    auto& singleton_proxy =
        Parallel::get_parallel_component<MemoryMonitor<Metavariables>>(cache);

    auto& group_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const int my_proc = Parallel::my_proc(*Parallel::local_branch(group_proxy));

    const double size_in_bytes = static_cast<double>(
        size_of_object_in_bytes(*Parallel::local_branch(group_proxy)));

    const double size_in_MB = size_in_bytes / 1.0e6;

    // See note in ProcessNodegroup on why we don't call
    // Parallel::contribute_to_reduction here
    Parallel::simple_action<ContributeMemoryData<ParallelComponent>>(
        singleton_proxy, time, my_proc, size_in_MB);
  }
};

/*!
 * \brief Simple action meant to be run on a Singleton that writes the size of
 * the Singleton component to disk.
 *
 * \details The columns in the dat file are
 *
 * - %Time
 * - Proc
 * - Size (MB)
 *
 * The dat file will be placed in the `/MemoryMonitors/` group in the reduction
 * file. The name of the dat file is the `pretty_type::name` of the component.
 */
struct ProcessSingleton {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/, const double time) {
    auto& singleton_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const double size_in_bytes = static_cast<double>(
        size_of_object_in_bytes(*Parallel::local(singleton_proxy)));
    const double size_in_MB = size_in_bytes / 1.0e6;

    const std::vector<std::string> legend{{"Time", "Proc", "Size (MB)"}};

    auto& observer_writer_proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);

    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        // Node 0 is always the writer
        observer_writer_proxy[0], subfile_name<ParallelComponent>(), legend,
        std::make_tuple(time,
                        Parallel::my_proc(*Parallel::local(singleton_proxy)),
                        size_in_MB));
  }
};
}  // namespace mem_monitor
