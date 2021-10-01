// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/OptionTags.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Time/Tags.hpp"
#include "Time/Time.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Tags {
/// \cond
struct ControlError;
struct Time;
}  // namespace Tags

namespace domain::OptionTags {
template <size_t>
struct DomainCreator;
}  // namespace domain::OptionTags
/// \endcond

namespace control_system {
namespace OptionHolder {
template <size_t DerivOrder>
struct ControlSystem;
}  // namespace OptionHolder

namespace OptionTags {
/// \ingroup OptionGroupsGroup
/// \ingroup ControlSystemGroup
struct ControlSystemGroup {
  static std::string name() noexcept { return "ControlSystems"; }
  static constexpr Options::String help = {"Control system"};
};

/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
template <typename ControlSystem>
struct ControlSystemInputs {
  static constexpr size_t deriv_order = ControlSystem::deriv_order;
  using type = OptionHolder::ControlSystem<deriv_order>;
  static constexpr Options::String help{"Options for a control system."};
  static std::string name() noexcept { return ControlSystem::name(); }
  using group = ControlSystemGroup;
};
}  // namespace OptionTags

namespace Tags {
// UsedForName is a control system struct (i.e. conforming to
// control_system::Protocols::ControlSystem)
template <typename ControlSystem>
struct ControlSystemInputs : db::SimpleTag {
  static constexpr size_t deriv_order = ControlSystem::deriv_order;
  using type = OptionHolder::ControlSystem<deriv_order>;
  using option_tags =
      tmpl::list<OptionTags::ControlSystemInputs<ControlSystem>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) noexcept {
    return option;
  }
};

struct ControlSystemName : db::SimpleTag {
  using type = std::string;
};

struct MeasureTranslationResult : db::SimpleTag {
  // using type = ::domain::Tags::MeshVelocity<2>::type;
  using type = std::vector<double>;
};

struct PrintOutput : db::SimpleTag {
  using type = bool;
};

template <size_t DerivOrder>
struct Averager : db::SimpleTag {
  using type = ::Averager<DerivOrder>;
};

struct TimescaleTuner : db::SimpleTag {
  using type = ::TimescaleTuner;
};

template <size_t DerivOrder>
struct Controller : db::SimpleTag {
  using type = ::Controller<DerivOrder>;
};

/// The measurement timescales associated with
/// domain::Tags::FunctionsOfTime.  Each function of time associated
/// with a control system has a corresponding set of timescales here,
/// represented as `PiecewisePolynomial<0>` with the same components
/// as the function itself.
struct MeasurementTimescales : db::SimpleTag {
  using type = std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;

  static constexpr bool pass_metavariables = true;

  template <typename Metavariables>
  using option_tags =
      tmpl::list<domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
                 ::OptionTags::InitialTimeStep>;

  template <typename Metavariables>
  static type create_from_options(
      const std::unique_ptr<::DomainCreator<Metavariables::volume_dim>>&
          domain_creator,
      const double initial_time_step) noexcept {
    std::unordered_map<std::string,
                       std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
        timescales;
    for (const auto& function_of_time : domain_creator->functions_of_time()) {
      if (function_of_time.second->time_bounds()[1] ==
          std::numeric_limits<double>::infinity()) {
        // This function of time is not controlled by a control
        // system.  It is an analytic function or similar.
        continue;
      }
      const double function_initial_time =
          function_of_time.second->time_bounds()[0];
      const DataVector used_for_size =
          function_of_time.second->func(function_initial_time)[0];

      // This check is intentionally inside the loop over the
      // functions of time so that it will not trigger for domains
      // without control systems.
      if (initial_time_step <= 0.0) {
        ERROR(
            "Control systems can only be used in forward-in-time evolutions.");
      }

      auto initial_timescale =
          make_with_value<DataVector>(used_for_size, initial_time_step);
      // For now just have constant measurements every time step
      timescales.emplace(
          function_of_time.first,
          std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
              function_initial_time, std::array{std::move(initial_timescale)},
              std::numeric_limits<double>::infinity()));
    }
    return timescales;
  }
};
}  // namespace Tags

namespace OptionHolder {
template <size_t DerivOrder>
struct ControlSystem {
  struct Averager {
    using type = ::Averager<DerivOrder>;
    static constexpr Options::String help = {"Options for the averager."};
  };

  struct Controller {
    using type = ::Controller<DerivOrder>;
    static constexpr Options::String help = {"Options for the controller."};
  };

  struct TimescaleTuner {
    using type = ::TimescaleTuner;
    static constexpr Options::String help = {
        "Options for the timescale tuner."};
  };

  struct PrintOutput {
    using type = bool;
    static constexpr Options::String help = {
        "Whether to print control system output or not."};
  };

  using options = tmpl::list<Averager, Controller,
                             TimescaleTuner, PrintOutput>;
  static constexpr Options::String help = {"Options for a controlsystem."};

  using db_tags = tmpl::list<control_system::Tags::Averager<DerivOrder>,
                             control_system::Tags::Controller<DerivOrder>,
                             control_system::Tags::TimescaleTuner,
                             control_system::Tags::PrintOutput>;
  using tuple_of_inputs = tuples::tagged_tuple_from_typelist<db_tags>;

  ControlSystem(::Averager<DerivOrder> averager,
                ::Controller<DerivOrder> controller, ::TimescaleTuner tuner,
                bool print_output) noexcept {
    tuples::get<control_system::Tags::Averager<DerivOrder>>(all_inputs_) =
        std::move(averager);
    tuples::get<control_system::Tags::Controller<DerivOrder>>(all_inputs_) =
        std::move(controller);
    tuples::get<control_system::Tags::TimescaleTuner>(all_inputs_) =
        std::move(tuner);
    tuples::get<control_system::Tags::PrintOutput>(all_inputs_) = print_output;
  }

  ControlSystem() = default;
  ControlSystem(const ControlSystem& /*rhs*/) = default;
  ControlSystem& operator=(const ControlSystem& /*rhs*/) = delete;
  ControlSystem(ControlSystem&& /*rhs*/) noexcept = default;
  ControlSystem& operator=(ControlSystem&& /*rhs*/) noexcept = default;
  ~ControlSystem() = default;

  // clang-tidy non-const reference pointer.
  void pup(PUP::er& p) noexcept { p | all_inputs_; };  // NOLINT

  tuple_of_inputs get_inputs() noexcept { return std::move(all_inputs_); }

 private:
  tuple_of_inputs all_inputs_{};
};

// template <size_t DerivOrder>
// bool operator==(const ControlSystem<DerivOrder>& lhs,
//                const ControlSystem<DerivOrder>& rhs) noexcept;
// template <size_t DerivOrder>
// bool operator!=(const ControlSystem<DerivOrder>& lhs,
//                const ControlSystem<DerivOrder>& rhs) noexcept;
}  // namespace OptionHolder
}  // namespace control_system
