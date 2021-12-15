// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>

#include "ControlSystem/InitialExpirationTimes.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/OptionTags.hpp"
#include "Time/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
template <class Metavariables, typename ControlSystem>
struct ControlComponent;
/// \endcond

namespace control_system::Tags {
/// \ingroup ControlSystemGroup
/// The FunctionsOfTime initialized from a DomainCreator, initial time
/// step, and control system OptionHolders.
struct FunctionsOfTimeInitialize : domain::Tags::FunctionsOfTime,
                                   db::SimpleTag {
  using type = std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;

  static constexpr bool pass_metavariables = true;

  static std::string name() { return "FunctionsOfTime"; }

  template <typename Metavariables>
  using option_tags = tmpl::flatten<
      tmpl::list<domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
                 ::OptionTags::InitialTime, ::OptionTags::InitialTimeStep,
                 control_system::inputs<tmpl::transform<
                     tmpl::filter<typename Metavariables::component_list,
                                  tt::is_a_lambda<ControlComponent, tmpl::_1>>,
                     tmpl::bind<tmpl::back, tmpl::_1>>>>>;

  template <typename Metavariables, typename... OptionHolders>
  static type create_from_options(
      const std::unique_ptr<::DomainCreator<Metavariables::volume_dim>>&
          domain_creator,
      const double initial_time, const double initial_time_step,
      const OptionHolders&... option_holders) {
    const std::unordered_map<std::string, double> initial_expiration_times =
        control_system::initial_expiration_times(
            initial_time, initial_time_step, option_holders...);

    auto functions_of_time =
        domain_creator->functions_of_time(initial_expiration_times);

    // Check that all control systems are actually controlling a function of
    // time, and that the expiration times have been set appropriately. If there
    // exists a control system that isn't controlling a function of time, or the
    // expiration times were set improperly, this is an error and we shouldn't
    // continue.
    for (const auto& [name, expr_time] : initial_expiration_times) {
      if (functions_of_time.count(name) == 0) {
        ERROR(
            "The control system '"
            << name
            << "' is not controlling a function of time. Check that the "
               "DomainCreator you have chosen uses all of the control "
               "systems in the executable. The existing functions of time are: "
            << keys_of(functions_of_time));
      }

      if (functions_of_time.at(name)->time_bounds()[1] != expr_time) {
        ERROR("The expiration time for the function of time '"
              << name << "' has been set improperly. It is supposed to be "
              << expr_time << " but is currently set to "
              << functions_of_time.at(name)->time_bounds()[1]
              << ". It is possible that the DomainCreator you are using isn't "
                 "compatible with the control systems.");
      }
    }

    return functions_of_time;
  }
};

namespace detail {

CREATE_HAS_STATIC_MEMBER_VARIABLE(override_functions_of_time)
CREATE_HAS_STATIC_MEMBER_VARIABLE_V(override_functions_of_time)

template <typename Metavariables, bool HasOverrideCubicFunctionsOfTime>
struct OptionList {
  using type = tmpl::conditional_t<
      Metavariables::override_functions_of_time,
      tmpl::flatten<tmpl::list<
          domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
          domain::FunctionsOfTime::OptionTags::FunctionOfTimeFile,
          domain::FunctionsOfTime::OptionTags::FunctionOfTimeNameMap,
          ::OptionTags::InitialTime, ::OptionTags::InitialTimeStep,
          control_system::inputs<tmpl::transform<
              tmpl::filter<typename Metavariables::component_list,
                           tt::is_a_lambda<ControlComponent, tmpl::_1>>,
              tmpl::bind<tmpl::back, tmpl::_1>>>>>,
      tmpl::list<domain::OptionTags::DomainCreator<Metavariables::volume_dim>>>;
};

template <typename Metavariables>
struct OptionList<Metavariables, false> {
  using type =
      tmpl::list<domain::OptionTags::DomainCreator<Metavariables::volume_dim>>;
};
}  // namespace detail

/// Tag to retreive the FunctionsOfTime from the GlobalCache.
struct FunctionsOfTime : db::BaseTag {};

/// \brief The FunctionsOfTime initialized from a DomainCreator or
/// (if `override_functions_of_time` is true in the metavariables) read
/// from a file.
///
/// \details When `override_functions_of_time == true` in the
/// metavariables, after obtaining the FunctionsOfTime from the DomainCreator,
/// one or more of those FunctionsOfTime (which must be cubic piecewise
/// polynomials) is overriden using data read from an HDF5 file via
/// domain::Tags::read_spec_piecewise_polynomial()
struct FunctionsOfTimeSpeCAndSpECTREControlSys : domain::Tags::FunctionsOfTime,
                                                 db::SimpleTag {
  using type = std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;

  static constexpr bool pass_metavariables = true;

  static std::string name() { return "FunctionsOfTime"; }

  template <typename Metavariables>
  using option_tags = typename detail::OptionList<
      Metavariables,
      ::detail::has_override_functions_of_time_v<Metavariables>>::type;

  template <typename Metavariables, typename... OptionHolders>
  static type create_from_options(
      const std::unique_ptr<::DomainCreator<Metavariables::volume_dim>>&
          domain_creator,
      const std::optional<std::string>& function_of_time_file,
      const std::optional<std::map<std::string, std::string>>&
          function_of_time_name_map,
      const double initial_time, const double initial_time_step,
      const OptionHolders&... option_holders) {
    if (function_of_time_file and function_of_time_name_map) {
      const auto initial_expiration_times =
          control_system::initial_expiration_times(
              initial_time, initial_time_step, option_holders...);

      auto functions_of_time =
          domain_creator->functions_of_time(initial_expiration_times);

      // Currently, only support order 2 or 3 piecewise polynomials.
      // This could be generalized later, but the SpEC functions of time
      // that we will read in with this action will always be 2nd-order or
      // 3rd-order piecewise polynomials
      std::unordered_map<std::string,
                         domain::FunctionsOfTime::PiecewisePolynomial<2>>
          spec_functions_of_time_second_order{};
      std::unordered_map<std::string,
                         domain::FunctionsOfTime::PiecewisePolynomial<3>>
          spec_functions_of_time_third_order{};

      // Import those functions of time of each supported order
      domain::FunctionsOfTime::read_spec_piecewise_polynomial(
          make_not_null(&spec_functions_of_time_second_order),
          *function_of_time_file, *function_of_time_name_map);
      domain::FunctionsOfTime::read_spec_piecewise_polynomial(
          make_not_null(&spec_functions_of_time_third_order),
          *function_of_time_file, *function_of_time_name_map);

      for (const auto& [spec_name, spectre_name] : *function_of_time_name_map) {
        (void)spec_name;
        // The FunctionsOfTime we are mutating must already have
        // an element with key==spectre_name; this action only
        // mutates the value associated with that key
        if (functions_of_time.count(spectre_name) == 0) {
          ERROR("Trying to import data for key "
                << spectre_name
                << " in FunctionsOfTime, but FunctionsOfTime does not "
                   "contain that key. This might happen if the option "
                   "FunctionOfTimeNameMap is not specified correctly. Keys "
                   "contained in FunctionsOfTime: "
                << keys_of(functions_of_time) << "\n");
        }
        auto* piecewise_polynomial_second_order =
            dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<2>*>(
                functions_of_time[spectre_name].get());
        auto* piecewise_polynomial_third_order =
            dynamic_cast<domain::FunctionsOfTime::PiecewisePolynomial<3>*>(
                functions_of_time[spectre_name].get());
        if (piecewise_polynomial_second_order == nullptr) {
          if (piecewise_polynomial_third_order == nullptr) {
            ERROR(
                "The function of time with name "
                << spectre_name
                << " is not a PiecewisePolynomial<2> or PiecewisePolynomial<3> "
                   "and so cannot be set using "
                   "read_spec_piecewise_polynomial\n");
          } else {
            *piecewise_polynomial_third_order =
                spec_functions_of_time_third_order.at(spectre_name);
          }
        } else {
          *piecewise_polynomial_second_order =
              spec_functions_of_time_second_order.at(spectre_name);
        }
      }

      return functions_of_time;
    } else {
      return domain_creator->functions_of_time();
    }
  }

  template <typename Metavariables>
  static type create_from_options(
      const std::unique_ptr<::DomainCreator<Metavariables::volume_dim>>&
          domain_creator) {
    return domain_creator->functions_of_time();
  }
};
}  // namespace control_system::Tags
