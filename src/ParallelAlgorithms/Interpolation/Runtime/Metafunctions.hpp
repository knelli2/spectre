// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/TagTraits.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp2::metafunctions {
namespace detail {
template <typename ComputeTag>
struct get_simple_tag {
  using type = typename ComputeTag::base;
};

template <typename ListOfTags>
struct simple_tags_from_mixed_tags {
  using partition = tmpl::partition<ListOfTags, db::is_compute_tag<tmpl::_1>>;
  using compute_to_simple_tags =
      tmpl::transform<typename partition::first_type, get_simple_tag<tmpl::_1>>;
  using type =
      tmpl::append<typename partition::second_type, compute_to_simple_tags>;
};
}  // namespace detail

template <typename ListOfTags>
using simple_tags_from_mixed_tags =
    typename detail::simple_tags_from_mixed_tags<ListOfTags>::type;

namespace detail {
template <typename Callback>
struct get_tags_to_observe {
  using type = typename Callback::tags_to_observe_on_target;
};
}  // namespace detail

template <typename Target>
using all_callbacks = tmpl::append<typename Target::possible_runtime_callbacks,
                                   typename Target::compile_time_callbacks>;

template <typename Target>
using all_tags_to_observe =
    simple_tags_from_mixed_tags<tmpl::flatten<tmpl::transform<
        all_callbacks<Target>, detail::get_tags_to_observe<tmpl::_1>>>>;
}  // namespace intrp2::metafunctions
