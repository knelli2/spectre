// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependentOptions/ExpansionMap.hpp"

#include <array>
#include <string>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace domain::creators::time_dependent_options {
ExpansionMapOptions::ExpansionMapOptions(
    const std::variant<std::array<double, 3>, FromVolumeFile<names::Expansion>>&
        expansion_values,
    const std::variant<std::array<double, 3>, FromVolumeFile<names::Expansion>>&
        expansion_outer_boundary_values,
    const double decay_timescale_outer_boundary_in,
    std::optional<double> decay_timescale_in,
    std::optional<double> asymptotic_velocity_outer_boundary_in,
    const Options::Context& context)
    : decay_timescale_outer_boundary(decay_timescale_outer_boundary_in),
      decay_timescale(std::move(decay_timescale_in)),
      asymptotic_velocity_outer_boundary(
          std::move(asymptotic_velocity_outer_boundary_in)) {
  if (decay_timescale.has_value() ==
      asymptotic_velocity_outer_boundary.has_value()) {
    PARSE_ERROR(context,
                "Must specify one of DecayTimescale or "
                "AsymptoticVelocityOuterBoundary, but not both.");
  }

  const auto set_values =
      [](const gsl::not_null<std::array<DataVector, 3>*> to_set,
         const auto& input_values, const bool is_outer_boundary) {
        if (std::holds_alternative<std::array<double, 3>>(input_values)) {
          auto& values = std::get<std::array<double, 3>>(input_values);
          for (size_t i = 0; i < to_set->size(); i++) {
            gsl::at(*to_set, i) = DataVector{1, gsl::at(values, i)};
          }
        } else if (std::holds_alternative<FromVolumeFile<names::Expansion>>(
                       input_values)) {
          auto& values_from_file =
              std::get<FromVolumeFile<names::Expansion>>(input_values);
          *to_set = is_outer_boundary
                        ? values_from_file.expansion_values_outer_boundary
                        : values_from_file.expansion_values;
        }
      };

  // Expansion values
  set_values(make_not_null(&initial_values), expansion_values, false);
  // Outer boundary
  set_values(make_not_null(&initial_values_outer_boundary),
             expansion_outer_boundary_values, true);
}
}  // namespace domain::creators::time_dependent_options
