// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <ostream>

#include "DataStructures/DataBox/Tag.hpp"

/// \cond
namespace Options {
class Option;
template <typename T>
struct create_from_yaml;
}  // namespace Options
/// \endcond

/*!
 * \brief State of an element
 */
enum class ElementState : int {
  ChuggingAlong = 0,
  WaitingForBoundaryDataInAction = 1,
  WaitingForBoundaryDataInPostProcessor = 2,
  WaitingForFoTInAction = 3,
  WaitingForFoTInPostProcessor = 4,
  WaitingForMeasurementTimescales = 5,
};
std::ostream& operator<<(std::ostream& os, ElementState ordering);

template <>
struct Options::create_from_yaml<ElementState> {
  template <typename Metavariables>
  static ElementState create(const Options::Option& options) {
    return create<void>(options);
  }
};

template <>
ElementState Options::create_from_yaml<ElementState>::create<void>(
    const Options::Option& options);

namespace Tags {
struct ElementState : db::SimpleTag {
  using type = ::ElementState;
};
}  // namespace Tags
