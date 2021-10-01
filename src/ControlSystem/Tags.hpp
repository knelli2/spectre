// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
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

namespace control_system::Tags {
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
                 OptionTags::InitialTimeStep>;

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

struct ControlSystemName : db::SimpleTag {
  using type = std::string;
};
}  // namespace control_system::Tags

namespace control_system::OptionHolder {
/// Options for control systems.
template <typename DerivOrder>
struct ControlSystem {
  struct Averager {
    using type = ::Averager<DerivOrder>;
    static constexpr Options::String help = {"Averager"};
    using group = OptionTags::ControlSystemGroup;
  };

  struct Controller {
    using type = ::Controller<DerivOrder>;
    static constexpr Options::String help = {"Controller"};
    using group = OptionTags::ControlSystemGroup;
  };

  struct TimescaleTuner {
    using type = ::TimescaleTuner;
    static constexpr Options::String help = {"TimescaleTuner"};
    using group = OptionTags::ControlSystemGroup;
  };

  struct PrintOutput {
    using type = bool;
    static constexpr Options::String help = {"PrintOutput"};
    using group = OptionTags::ControlSystemGroup;
  };

  using options = tmpl::list<Averager<DerivOrder>, Controller<DerivOrder>,
                             TimescaleTuner, PrintOuput>;
  static constexpr Options::String help = {"Options for a controlsystem."};

  ControlSystem(Averager<DerivOrder> averager,
                Controller<DerivOrder> controller, TimescaleTuner tuner,
                bool PrintOutput) noexcept;

  ControlSystem() = default;
  ControlSystem(const ApparentHorizon& /*rhs*/) = default;
  ControlSystem& operator=(const ApparentHorizon& /*rhs*/) = delete;
  ControlSystem(ApparentHorizon&& /*rhs*/) noexcept = default;
  ControlSystem& operator=(ApparentHorizon&& /*rhs*/) noexcept = default;
  ~ControlSystem() = default;

  // clang-tidy non-const reference pointer.
  void pup(PUP::er& p) noexcept;  // NOLINT

  Averager<DerivOrder> averager_;
  Controller<DerivOrder> controller_;
  Timescaletuner tuner_;
  bool print_output_;
};

template <typename DerivOrder>
bool operator==(const ControlSystem<DerivOrder>& lhs,
                const ControlSystem<DerivOrder>& rhs) noexcept;
template <typename DerivOrder>
bool operator!=(const ControlSystem<DerivOrder>& lhs,
                const ControlSystem<DerivOrder>& rhs) noexcept;

}  // namespace OptionHolders

namespace OptionTags {

struct ControlSystemGroup {
  static std::string name() noexcept { return "ControlSystem"; }
  static constexpr Options::String help = {"Control system"};
};
}  // namespace OptionTags

template <typename DerivOrder>
struct ControlSystem {
  using type = OptionHolders::ControlSystem<DerivOrder>;
  static constexpr Options::String help{
      "Options for interpolation onto apparent horizon."};
  static std::string name() noexcept {
    return Options::name<InterpolationTargetTag>();
  }
  using group = ControlSystemGroup;
};

namespace Tags {
template <typename DerivOrder>
struct ControlSystem : db::SimpleTag {
  using type = OptionHolders::ControlSystem<DerivOrder>;
  using option_tags =
      tmpl::list<OptionTags::ControlSystem<InterpolationTargetTag, DerivOrder>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) noexcept {
    return option;
  }
};
}

// This is alpha_m in Dan's paper
// struct MeasurementTimeScaleOverExpirationDeltaT {
//  using type = double;
//  static constexpr Options::String help = {
//      "How often to measure control error as a fraction of the expiration "
//      "timescale of the function of time"};
//  using group = OptionTags::ControlSystemGroup;
//};
//
//// This is alpha_d in Dan's paper
// struct ExpirationDeltaTOverDampingTimescale {
//  using type = double;
//  static constexpr Options::String help = {
//      "How long function of time expiration should be as a fraction of the "
//      "damping timescale"};
//  using group = OptionTags::ControlSystemGroup;
//};

namespace Tags {

// struct MeasurementTimeScaleOverExpirationDeltaT : db::SimpleTag {
//  using type = double;
//  using option_tags =
//      tmpl::list<OptionTags::MeasurementTimeScaleOverExpirationDeltaT>;
//  static constexpr bool pass_metavariables = false;
//  static auto create_from_options(const double value) { return value; }
//};
//
// struct ExpirationDeltaTOverDampingTimescale : db::SimpleTag {
//  using type = double;
//  using option_tags =
//      tmpl::list<OptionTags::ExpirationDeltaTOverDampingTimescale>;
//  static constexpr bool pass_metavariables = false;
//  static auto create_from_options(const double value) { return value; }
//};

struct MeasureTranslationResult : db::SimpleTag {
  // using type = ::domain::Tags::MeshVelocity<2>::type;
  using type = std::vector<double>;
};

struct PrintOutput : db::SimpleTag {
  using type = bool;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::PrintOutput>;

  static auto create_from_options(const bool printoutput) noexcept {
    return printoutput;
  }
};

template <size_t DerivOrder>
struct Averager : db::SimpleTag {
  using type = ::Averager<DerivOrder>;

  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::Averager<DerivOrder>>;

  static auto create_from_options(
      const ::Averager<DerivOrder>& averager) noexcept {
    return averager;
  }
};

struct TimescaleTuner : db::SimpleTag {
  using type = ::TimescaleTuner;

  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::TimescaleTuner>;

  static auto create_from_options(
      const ::TimescaleTuner& timescale_tuner) noexcept {
    return timescale_tuner;
  }
};

template <size_t DerivOrder>
struct Controller : db::SimpleTag {
  using type = ::Controller<DerivOrder>;

  static constexpr bool pass_metavariables = false;

  using option_tags = tmpl::list<OptionTags::Controller<DerivOrder>>;

  static auto create_from_options(
      const ::Controller<DerivOrder>& controller) noexcept {
    return controller;
  }
};

}  // namespace Tags
