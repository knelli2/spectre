// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <ostream>

/// \cond
namespace Options {
class Option;
template <typename T>
struct create_from_yaml;
}  // namespace Options
/// \endcond

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {

/*!
 * \brief Type of falloff that can be used to determine which
 * `ShapeMapTransitionFunction` to choose.
 *
 * \details Main use is for options in the input file
 */
enum class Falloff { Linear, Inverse };

std::ostream& operator<<(std::ostream& os, Falloff falloff);

}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions

template <>
struct Options::create_from_yaml<
    domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff> {
  template <typename Metavariables>
  static domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff create(
      const Options::Option& options) {
    return create<void>(options);
  }
};
template <>
domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff
Options::create_from_yaml<
    domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff>::
    create<void>(const Options::Option& options);
