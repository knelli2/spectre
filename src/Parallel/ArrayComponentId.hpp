// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/functional/hash.hpp>
#include <ckarrayindex.h>
#include <cstddef>
#include <functional>
#include <string>
#include <type_traits>

#include "Parallel/ArrayIndex.hpp"
#include "Parallel/StringIndex.hpp"
#include "Utilities/PrettyType.hpp"

namespace PUP {
class er;
}  // namespace PUP

namespace Parallel {
/*!
 * \brief An ID type that identifies both the parallel component and the index
 * in the parallel component.
 *
 * A specialization of `std::hash` is provided to allow using `ArrayComponentId`
 * as a key in associative containers.
 */
class ArrayComponentId {
 public:
  ArrayComponentId() = default;  // Needed for Charm++ serialization

  template <typename ParallelComponent>
  ArrayComponentId(const ParallelComponent* /*meta*/,
                   const CkArrayIndex& index);

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  size_t component_id() const { return component_id_; }

  const CkArrayIndex& array_index() const { return array_index_; }

 private:
  size_t component_id_{0};
  CkArrayIndex array_index_{};
};

template <typename ParallelComponent>
ArrayComponentId::ArrayComponentId(const ParallelComponent* const /*meta*/,
                                   const CkArrayIndex& index)
    : component_id_(
          std::hash<std::string>{}(pretty_type::get_name<ParallelComponent>())),
      array_index_(index) {}

bool operator==(const ArrayComponentId& lhs, const ArrayComponentId& rhs);

bool operator!=(const ArrayComponentId& lhs, const ArrayComponentId& rhs);

std::ostream& operator<<(std::ostream& os,
                         const ArrayComponentId& array_component_id);

/*!
 * \brief A convenience function that will make an `ArrayComponentId` from the
 * templated `ParallelComponent` and the passed in `array_index`.
 */
template <typename ParallelComponent, typename ArrayIndexType>
ArrayComponentId make_array_component_id(const ArrayIndexType& array_index) {
  return ArrayComponentId{std::add_pointer_t<ParallelComponent>{nullptr},
                          Parallel::ArrayIndex<ArrayIndexType>(array_index)};
}

template <typename ParallelComponent>
ArrayComponentId make_array_component_id(const std::string& array_index) {
  return ArrayComponentId{
      std::add_pointer_t<ParallelComponent>{nullptr},
      Parallel::ArrayIndex<StringIndex>(StringIndex{array_index})};
}
}  // namespace Parallel

namespace std {
template <>
struct hash<Parallel::ArrayComponentId> {
  size_t operator()(const Parallel::ArrayComponentId& t) const {
    size_t result = t.component_id();
    boost::hash_combine(result, t.array_index().hash());
    return result;
  }
};
}  // namespace std
