// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/OptionTags.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace domain::Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup ComputationalDomainGroup
/// Base tag to retrieve the centers of objects in the domain corresponding to
/// the `ObjectLabel`.
template <ObjectLabel Label>
struct ObjectCenter : db::BaseTag {};

/// \ingroup DataBoxTagsGroup
/// \ingroup ComputationalDomainGroup
/// The center of the excision sphere for the given object in the given frame.
///
/// Even though this can easily be retrieved from the domain, we add it as its
/// own tag so we can access it through the base tag. This way, other things
/// (like the control system) can grab the center and be agnostic to what the
/// object actually is.
template <ObjectLabel Label, typename Frame>
struct ExcisionCenter : ObjectCenter<Label>, db::SimpleTag {
  using type = tnsr::I<double, 3, Frame>;
  static std::string name() {
    return pretty_type::name<Frame>() + "CenterObject" + get_output(Label);
  }

  using option_tags = tmpl::list<domain::OptionTags::DomainCreator<3>>;
  static constexpr bool pass_metavariables = false;

  static type create_from_options(
      const std::unique_ptr<::DomainCreator<3>>& domain_creator);
};
}  // namespace domain::Tags
