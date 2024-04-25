// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>
#include <string>

class StringIndex {
 public:
  StringIndex() = default;
  StringIndex(const StringIndex&) = default;
  StringIndex& operator=(const StringIndex&) = default;
  StringIndex(StringIndex&&) = default;
  StringIndex& operator=(StringIndex&&) = default;
  ~StringIndex() = default;

  StringIndex(const std::string& name);

  // NOLINTNEXTLINE
  friend bool operator==(const StringIndex& index_1,
                         const StringIndex& index_2);

 private:
  std::uint32_t hash_of_name_;
  // #pragma GCC diagnostic push
  // #pragma GCC diagnostic ignored "-Wunused-private-field"
  std::uint32_t padding_1_;
  std::uint32_t padding_2_;
  // #pragma GCC diagnostic pop
  static_assert(sizeof(hash_of_name_) + sizeof(padding_1_) +
                    sizeof(padding_2_) ==
                3 * sizeof(int));
};

bool operator!=(const StringIndex& index_1, const StringIndex& index_2);

// clang-tidy: do not modify namespace std
namespace std {  // NOLINT
template <>
struct hash<StringIndex> {
  size_t operator()(const StringIndex& index) const;
};
}  // namespace std

/// \cond
PUPbytes(StringIndex)  // NOLINT
                       /// \endcond
