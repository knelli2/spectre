// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <utility>

#include "ApparentHorizons/ObjectLabel.hpp"
#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "Options/Options.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace control_system {
/// \cond
template <typename ControlSystem>
struct OptionHolder;

namespace ah {
enum class HorizonLabel;
}
/// \endcond

namespace QueueTags {
template <::ah::ObjectLabel WhichHorizon>
struct Center {
  using type = DataVector;
};

template <typename Frame>
struct Strahlkorper {
  using type = ::Strahlkorper<Frame>;
};
}  // namespace QueueTags

/// \ingroup ControlSystemGroup
/// All option tags related to the control system
namespace OptionTags {
/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
/// Options group for all control system options
struct ControlSystemGroup {
  static std::string name() { return "ControlSystems"; }
  static constexpr Options::String help = {
      "Options for all control systems used in a simulation."};
};

/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
/// Option tag for each individual control system. The name of this option is
/// the name of the \p ControlSystem struct it is templated on. This way all
/// control systems will have a unique name.
template <typename ControlSystem>
struct ControlSystemInputs {
  using type = control_system::OptionHolder<ControlSystem>;
  static constexpr Options::String help{"Options for a control system."};
  static std::string name() { return ControlSystem::name(); }
  using group = ControlSystemGroup;
};

/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
/// Option tag on whether to write data to disk.
struct WriteDataToDisk {
  using type = bool;
  static constexpr Options::String help = {
      "Whether control system data should be saved during an evolution."};
  using group = ControlSystemGroup;
};

struct RestrictToRotationAboutZAxis {
  using type = bool;
  static constexpr Options::String help = {
      "Restrict rotation to be about the Z axis."};
  using group = ControlSystemGroup;
};
}  // namespace OptionTags

/// \ingroup ControlSystemGroup
/// Alias to get all the option holders from a list of control systems. This is
/// useful in the `option_tags` alias of simple tags for getting all the options
/// from control systems.
template <typename ControlSystems>
using inputs =
    tmpl::transform<ControlSystems,
                    tmpl::bind<OptionTags::ControlSystemInputs, tmpl::_1>>;

/// \ingroup ControlSystemGroup
/// All DataBox tags related to the control system
namespace Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the name of a control system
struct ControlSystemName : db::SimpleTag {
  using type = std::string;
};

struct RestrictToRotationAboutZAxis : db::SimpleTag {
  using type = bool;
  using option_tags = tmpl::list<OptionTags::RestrictToRotationAboutZAxis>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for writing control system data to disk
struct WriteDataToDisk : db::SimpleTag {
  using type = bool;
  using option_tags = tmpl::list<OptionTags::WriteDataToDisk>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for all options of a single control system.
///
/// Only intended to be used during the initialization phase as a way of getting
/// options from multiple control systems into their corresponding components
/// DataBox.
template <typename ControlSystem>
struct ControlSystemInputs : db::SimpleTag {
  using type = control_system::OptionHolder<ControlSystem>;
  using option_tags =
      tmpl::list<OptionTags::ControlSystemInputs<ControlSystem>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the averager
template <size_t DerivOrder>
struct Averager : db::SimpleTag {
  using type = ::Averager<DerivOrder>;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the timescale tuner
struct TimescaleTuner : db::SimpleTag {
  using type = ::TimescaleTuner;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the controller
template <size_t DerivOrder>
struct Controller : db::SimpleTag {
  using type = ::Controller<DerivOrder>;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the control error
template <typename ControlSystem>
struct ControlError : db::SimpleTag {
  using type = typename ControlSystem::control_error;
};
}  // namespace Tags

/// \ingroup ControlSystemGroup
/// Holds all options for a single control system
///
/// This struct collects all the options for a given control system during
/// option parsing. Then during initialization, the options can be retrieved via
/// their public member names and assigned to their corresponding DataBox tags.
template <typename ControlSystem>
struct OptionHolder {
  static_assert(tt::assert_conforms_to<
                ControlSystem, control_system::protocols::ControlSystem>);
  using control_system = ControlSystem;
  static constexpr size_t deriv_order = control_system::deriv_order;
  struct Averager {
    using type = ::Averager<deriv_order - 1>;
    static constexpr Options::String help = {
        "Averages the derivatives of the control error and possibly the "
        "control error itself."};
  };

  struct Controller {
    using type = ::Controller<deriv_order>;
    static constexpr Options::String help = {
        "Computes the control signal which will be used to reset the functions "
        "of time."};
  };

  struct TimescaleTuner {
    using type = ::TimescaleTuner;
    static constexpr Options::String help = {
        "Keeps track of the damping timescales for the control system upon "
        "which other timescales are based of off."};
  };

  struct ControlError {
    using type = typename ControlSystem::control_error;
    static constexpr Options::String help = {
        "Computes the control error for the control system based on quantities "
        "in the simulation."};
  };

  using options =
      tmpl::list<Averager, Controller, TimescaleTuner, ControlError>;
  static constexpr Options::String help = {"Options for a control system."};

  OptionHolder(::Averager<deriv_order - 1> input_averager,
               ::Controller<deriv_order> input_controller,
               ::TimescaleTuner input_tuner,
               typename ControlSystem::control_error input_control_error)
      : averager(std::move(input_averager)),
        controller(std::move(input_controller)),
        tuner(std::move(input_tuner)),
        control_error(std::move(input_control_error)) {}

  OptionHolder() = default;
  OptionHolder(const OptionHolder& /*rhs*/) = default;
  OptionHolder& operator=(const OptionHolder& /*rhs*/) = default;
  OptionHolder(OptionHolder&& /*rhs*/) = default;
  OptionHolder& operator=(OptionHolder&& /*rhs*/) = default;
  ~OptionHolder() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | averager;
    p | controller;
    p | tuner;
    p | control_error;
  };

  // These members are specifically made public for easy access during
  // initialization
  ::Averager<deriv_order - 1> averager{};
  ::Controller<deriv_order> controller{};
  ::TimescaleTuner tuner{};
  typename ControlSystem::control_error control_error{};
};
}  // namespace control_system
