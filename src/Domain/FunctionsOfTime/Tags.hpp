// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>

#include "ControlSystem/Component.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/OptionTags.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/ReadSpecPiecewisePolynomial.hpp"
#include "Domain/OptionTags.hpp"
#include "Options/Options.hpp"
#include "Time/Tags.hpp"
#include "Utilities/StdHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/CreateHasStaticMemberVariable.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// \cond
template <size_t VolumeDim>
class DomainCreator;
/// \endcond

namespace detail {

CREATE_HAS_STATIC_MEMBER_VARIABLE(override_functions_of_time)
CREATE_HAS_STATIC_MEMBER_VARIABLE_V(override_functions_of_time)

template <typename Metavariables, bool HasOverrideCubicFunctionsOfTime>
struct OptionList {
  using type = tmpl::conditional_t<
      Metavariables::override_functions_of_time,
      tmpl::list<domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
                 domain::FunctionsOfTime::OptionTags::FunctionOfTimeFile,
                 domain::FunctionsOfTime::OptionTags::FunctionOfTimeNameMap>,
      tmpl::flatten<tmpl::list<
          domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
          ::OptionTags::InitialSlabSize,
          control_system::option_holders<tmpl::transform<
              tmpl::filter<typename Metavariables::component_list,
                           tt::is_a_lambda<ControlComponent, tmpl::_1>>,
              tmpl::bind<tmpl::back, tmpl::_1>>>>>>;
};

template <typename Metavariables>
struct OptionList<Metavariables, false> {
  using type = tmpl::flatten<
      tmpl::list<domain::OptionTags::DomainCreator<Metavariables::volume_dim>,
                 ::OptionTags::InitialSlabSize,
                 control_system::option_holders<tmpl::transform<
                     tmpl::filter<typename Metavariables::component_list,
                                  tt::is_a_lambda<ControlComponent, tmpl::_1>>,
                     tmpl::bind<tmpl::back, tmpl::_1>>>>>;
};
}  // namespace detail

namespace domain::Tags {
/// \brief The FunctionsOfTime initialized from a DomainCreator and control
/// systems or (if `override_functions_of_time` is true in the metavariables)
/// read from a file.
///
/// \details When `override_functions_of_time == true` in the
/// metavariables, after obtaining the FunctionsOfTime from the DomainCreator,
/// one or more of those FunctionsOfTime (which must be cubic piecewise
/// polynomials) is overriden using data read from an HDF5 file via
/// domain::Tags::read_spec_piecewise_polynomial()
///
/// When `override_functions_of_time == false`, some extra options are passed in
/// addition to the DomainCreator. First is the initial slab size. The
/// (potential) second and onward options are control system OptionHolders which
/// are used to set the initial expiration times. They are potential options
/// because the number of OptionHolders depends on the number of
/// ControlComponents in the `component_list` type alias of the metavariables.
struct FunctionsOfTime : db::SimpleTag {
  using type = std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;

  static constexpr bool pass_metavariables = true;

  template <typename Metavariables>
  using option_tags = typename ::detail::OptionList<
      Metavariables,
      ::detail::has_override_functions_of_time_v<Metavariables>>::type;

  template <typename Metavariables>
  static type create_from_options(
      const std::unique_ptr<::DomainCreator<Metavariables::volume_dim>>&
          domain_creator,
      const std::optional<std::string>& function_of_time_file,
      const std::optional<std::map<std::string, std::string>>&
          function_of_time_name_map) {
    if (function_of_time_file and function_of_time_name_map) {
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

      auto functions_of_time{domain_creator->functions_of_time()};
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

  template <typename Metavariables, typename... OptionHolders>
  static type create_from_options(
      const std::unique_ptr<::DomainCreator<Metavariables::volume_dim>>&
          domain_creator,
      const double initial_slab_size, const OptionHolders&&... option_holders) {
    std::unordered_map<std::string, double> initial_expiration_times{};
    [[maybe_unused]] const auto lambda = [&initial_expiration_times,
                                          &initial_slab_size](
                                             const auto& option_holder) {
      const auto controller = option_holder.controller;
      const std::string name = option_holder.name;
      const auto tuner = option_holder.tuner;

      const double update_fraction = controller.get_update_fraction();
      const double curr_timescale = min(tuner.current_timescale());
      const double initial_expiration_time = update_fraction * curr_timescale;
      initial_expiration_times[name] =
          initial_expiration_time < initial_slab_size ? initial_slab_size
                                                      : initial_expiration_time;
    };
    EXPAND_PACK_LEFT_TO_RIGHT(lambda(option_holders));
    return domain_creator->functions_of_time(initial_expiration_times);
  }
};
}  // namespace domain::Tags
