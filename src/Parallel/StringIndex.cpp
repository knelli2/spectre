// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Parallel/StringIndex.hpp"

#include <boost/functional/hash.hpp>
#include <string>

// The `static_assert` verifies that `StringIndex` satisfies the constraints
// imposed by Charm++ to make `StringIndex` able to act as an index into
// Charm++'s arrays. These constraints are:
// - `StringIndex` must satisfy `std::is_standard_layout` and `std::is_trivial`
// - `StringIndex` must not be larger than the size of three `int`s, i.e.
//   `sizeof(StringIndex) <= 3 * sizeof(int)`
static_assert(std::is_standard_layout_v<StringIndex> and
              std::is_trivial_v<StringIndex>);
static_assert(sizeof(StringIndex) == 3 * sizeof(int),
              "Wrong size for StringIndex");

StringIndex::StringIndex(const std::string& name)
    : hash_of_name_(std::hash<std::string>{}(name)),
      padding_1_{0},
      padding_2_{0} {}

bool operator==(const StringIndex& index_1, const StringIndex& index_2) {
  return index_1.hash_of_name_ == index_2.hash_of_name_;
}

bool operator!=(const StringIndex& index_1, const StringIndex& index_2) {
  return not(index_1 == index_2);
}

size_t hash_value(const StringIndex& index) {
  // StringIndex is used as an opaque array of bytes by Charm, so we
  // treat it that way as well.
  // clang-tidy: do not use reinterpret_cast
  const auto bytes = reinterpret_cast<const char*>(&index);  // NOLINT
  // clang-tidy: do not use pointer arithmetic
  return boost::hash_range(bytes, bytes + sizeof(index));  // NOLINT
}

// clang-tidy: do not modify namespace std
namespace std {  // NOLINT
size_t hash<StringIndex>::operator()(const StringIndex& index) const {
  return hash_value(index);
}
}  // namespace std
