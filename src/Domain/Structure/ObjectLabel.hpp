// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <iosfwd>
#include <string>

namespace domain {
/// Labels for the objects in a binary system.
enum class ObjectLabel : int {
  // HUGE WARNING IF YOU ARE GOING TO EDIT THIS ENUM!
  // Other parts of the code will expect the integral values for each of these
  // object labels to be exactly what they are here. If you feel you need to
  // change these, think very carefully about if you really need to. Consult a
  // core developer if you really feel you need to change the underlying values.
  // If it is decided that changing the values is the right thing to do, you
  // must find all other places in the code that assumed the previous values and
  // make changes there as well.

  /// The object along the positive x-axis in the grid frame
  A = 0,
  /// The object along the negative x-axis in the grid frame
  B = 1,
  /// A third object, typically centered at the origin.
  C = 2,
  /// This object has no label
  None = -1
};

std::string name(const ObjectLabel x);

std::ostream& operator<<(std::ostream& s, const ObjectLabel x);

/// \brief Similar to a `tmpl::list` but for `ObjectLabel`s.
template <ObjectLabel... Objects>
struct object_list {};
}  // namespace domain
