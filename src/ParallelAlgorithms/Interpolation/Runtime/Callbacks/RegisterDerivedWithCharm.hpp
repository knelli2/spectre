// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp2 {
namespace callbacks {
namespace detail {
template <typename Target>
struct get_callbacks {
  using type = tmpl::append<typename Target::possible_runtime_callbacks,
                            typename Target::compile_time_callbacks>;
};
}  // namespace detail

template <typename Targets>
void register_derived_with_charm() {
  using all_callbacks =
      tmpl::flatten<tmpl::transform<Targets, detail::get_callbacks<tmpl::_1>>>;
  register_classes_with_charm(all_callbacks{});
}
}  // namespace callbacks
}  // namespace intrp2
