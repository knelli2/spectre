// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/CreateIsCallable.hpp"

namespace Parallel {
namespace detail {
CREATE_IS_CALLABLE(output_inbox)
CREATE_IS_CALLABLE_V(output_inbox)
}  // namespace detail

template <typename InboxTag, typename... InboxTypes>
std::string output_inbox(const tuples::TaggedTuple<InboxTypes...>& inboxes,
                         const size_t indent_pad) {
  static_assert(tmpl::list_contains_v<tmpl::list<InboxTypes...>, InboxTag>);
  static_assert(detail::is_output_inbox_callable_v<
                    InboxTag, const typename InboxTag::type&, const size_t>,
                "To output an inbox, there must be a static 'output_inbox' "
                "function defined in the inbox tag.");

  return InboxTag::output_inbox(tuples::get<InboxTag>(inboxes), indent_pad);
}
}  // namespace Parallel
