// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <unordered_map>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Checkpoint.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags/SystemTags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/File.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/NodeLock.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system::Actions {
/*!
 * \brief Threaded action to be run on the 0th node of the ObserverWriter that
 * writes to the control system checkpoint subfile in the checkpoint H5 file.
 */
struct WriteSynchronousCheckpoint {
  template <typename ParallelComponent, typename DbTagsList,
            typename Metavariables, typename ArrayIndex, typename... Ts>
  static void apply(db::DataBox<DbTagsList>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const gsl::not_null<Parallel::NodeLock*> /*node_lock*/,
                    const std::string& subfile_name,
                    const std::vector<char>& averager,
                    const std::vector<char>& tuner,
                    const std::vector<char>& controller,
                    const std::vector<char>& control_error_class) {
    auto& reduction_file_lock =
        db::get_mutable_reference<observers::Tags::H5FileLock>(
            make_not_null(&box));
    const std::lock_guard hold_lock(reduction_file_lock);

    const std::string& h5_filename = cache.get_current_checkpoint_file();

    h5::H5File<h5::AccessType::ReadWrite> h5_file{h5_filename, true};

    auto& checkpoint_subfile =
        h5_file.try_insert<control_system::Checkpoint>(subfile_name);

    checkpoint_subfile.write_checkpoint(averager, tuner, controller,
                                        control_error_class);

    h5_file.close_current_object();
  }
};

/*!
 * \brief Iterable action that calls `WriteSynchronousCheckpoint` for the
 * \p ControlSystem which writes the control system checkpoint subfile.
 */
template <typename ControlSystem>
struct WriteCheckpoint {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    const auto& averager = db::get<Tags::Averager<ControlSystem>>(box);
    const auto& tuner = db::get<Tags::TimescaleTuner<ControlSystem>>(box);
    const auto& controller = db::get<Tags::Controller<ControlSystem>>(box);
    const auto& control_error_class =
        db::get<Tags::ControlError<ControlSystem>>(box);

    auto observer_proxy = Parallel::get_parallel_component<
        ::observers::ObserverWriter<Metavariables>>(cache)[0];

    Parallel::threaded_action<WriteSynchronousCheckpoint>(
        observer_proxy, "/ControlSystem/" + ControlSystem::name(),
        serialize(averager), serialize(tuner), serialize(controller),
        serialize(control_error_class));

    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace control_system::Actions
