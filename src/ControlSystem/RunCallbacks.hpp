// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Protocols/Submeasurement.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

#include "Parallel/Printf.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
/// \endcond

namespace control_system {
/// \ingroup ControlSystemGroup
/// Apply the `process_measurement` struct of each of the \p
/// ControlSystems to the result of the \p Submeasurement.
///
/// The submeasurement results are supplied as a `db::DataBox` in
/// order to allow individual control systems to select the portions
/// of the submeasurement that they are interested in.
///
/// In addition to being manually called, this struct is designed to
/// be usable as a `post_horizon_find_callback` from
/// `intrp::callbacks::FindApparentHorizon`.
template <typename Submeasurement, typename ControlSystems>
struct RunCallbacks {
 private:
  static_assert(
      tt::assert_conforms_to<Submeasurement, protocols::Submeasurement>);
  template <typename System>
  using assert_conforms_to_t = std::bool_constant<
      tt::assert_conforms_to<System, protocols::ControlSystem>>;
  static_assert(tmpl::all<ControlSystems,
                          tmpl::bind<assert_conforms_to_t, tmpl::_1>>::value);

 public:
  using observation_types = tmpl::list<>;

  template <typename DbTags, typename Metavariables>
  static void apply(const db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const LinkedMessageId<double>& measurement_id) noexcept {
    //Parallel::printf("Inside RunCallBacks\n");
    tmpl::for_each<ControlSystems>(
        [&box, &cache, &measurement_id](auto control_system_v) noexcept {
          using ControlSystem = tmpl::type_from<decltype(control_system_v)>;
          //Parallel::printf("  for_each control system inside RunCallBacks\n");
          db::apply<typename ControlSystem::process_measurement::
                        template argument_tags<Submeasurement>>(
              [&cache, &measurement_id](const auto&... args) noexcept {
                ControlSystem::process_measurement::apply(
                    Submeasurement{}, args..., cache, measurement_id);
              },
              box);
        });
  }
};
}  // namespace control_system
