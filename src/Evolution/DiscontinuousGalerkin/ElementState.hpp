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
  WaitingForMeasurementTimescales1 = -1,
  WaitingForMeasurementTimescales2 = -2,
  WaitingForMeasurementTimescales3 = -3,
  WaitingForMeasurementTimescales4 = -4,
  WaitingForMeasurementTimescales5 = -5,
  WaitingForMeasurementTimescales6 = -6,
  WaitingForMeasurementTimescales7 = -7,
  WaitingForMeasurementTimescales8 = -8,
  WaitingForMeasurementTimescales9 = -9,
  WaitingForMeasurementTimescales10 = -10,
  WaitingForMeasurementTimescales11 = -11,
  WaitingForMeasurementTimescales12 = -12,
  WaitingForMeasurementTimescales13 = -13,
  WaitingForMeasurementTimescales14 = -14,
  WaitingForMeasurementTimescales15 = -15,
  WaitingForMeasurementTimescales16 = -16,
  WaitingForMeasurementTimescales17 = -17,
  WaitingForMeasurementTimescales18 = -18,
  WaitingForMeasurementTimescales19 = -19,
  WaitingForMeasurementTimescales20 = -20,
  WaitingForMeasurementTimescales21 = -21,
  WaitingForMeasurementTimescales22 = -22,
  WaitingForMeasurementTimescales23 = -23,
  WaitingForMeasurementTimescales24 = -24,
  WaitingForMeasurementTimescales25 = -25,
  WaitingForMeasurementTimescales26 = -26,
  WaitingForMeasurementTimescales27 = -27,
  WaitingForMeasurementTimescales28 = -28,
  WaitingForMeasurementTimescales29 = -29,
  WaitingForMeasurementTimescales30 = -30,
  WaitingForMeasurementTimescales31 = -31,
  WaitingForMeasurementTimescales32 = -32,
  WaitingForMeasurementTimescales33 = -33,
  WaitingForMeasurementTimescales34 = -34,
  WaitingForMeasurementTimescales35 = -35,
  WaitingForMeasurementTimescales36 = -36,
  WaitingForMeasurementTimescales37 = -37,
  WaitingForMeasurementTimescales38 = -38,
  WaitingForMeasurementTimescales39 = -39,
  WaitingForMeasurementTimescales40 = -40,
  WaitingForMeasurementTimescales41 = -41,
  WaitingForMeasurementTimescales42 = -42,
  WaitingForMeasurementTimescales43 = -43,
  WaitingForMeasurementTimescales44 = -44,
  WaitingForMeasurementTimescales45 = -45,
  WaitingForMeasurementTimescales46 = -46,
  WaitingForMeasurementTimescales47 = -47,
  WaitingForMeasurementTimescales48 = -48,
  WaitingForMeasurementTimescales49 = -49,
  WaitingForMeasurementTimescales50 = -50,
  WaitingForMeasurementTimescales51 = -51,
  WaitingForMeasurementTimescales52 = -52,
  WaitingForMeasurementTimescales53 = -53,
  WaitingForMeasurementTimescales54 = -54,
  WaitingForMeasurementTimescales55 = -55,
  WaitingForMeasurementTimescales56 = -56,
  WaitingForMeasurementTimescales57 = -57,
  WaitingForMeasurementTimescales58 = -58,
  WaitingForMeasurementTimescales59 = -59,
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
