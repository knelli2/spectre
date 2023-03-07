// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/Tags/ObjectCenter.hpp"

#include <cstddef>
#include <memory>
#include <string>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/PrettyType.hpp"

namespace Frame {
struct Grid;
struct Distorted;
}  // namespace Frame

namespace domain::Tags {
template <ObjectLabel Label, typename Frame>
tnsr::I<double, 3, Frame> ExcisionCenter<Label, Frame>::create_from_options(
    const std::unique_ptr<::DomainCreator<3>>& domain_creator) {
  const auto domain = domain_creator->create_domain();
  // No frame here because the ExcisionSphere object is always in the Grid frame
  const std::string name = "ExcisionSphere"s + get_output(Label);
  if (domain.excision_spheres().count(name) != 1) {
    ERROR(name << " is not in the domains excision spheres but is needed to "
                  "generate the ExcisionCenter<"
               << Label << "," << pretty_type::name<Frame>() << ">.");
  }

  tnsr::I<double, 3, Frame> result{};
  tnsr::I<double, 3, ::Frame::Grid> center =
      domain.excision_spheres().at(name).center();

  for (size_t i = 0; i < 3; i++) {
    result.get(i) = center.get(i);
  }

  return result;
}

#define OBJECT(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) template struct ObjectCenter<OBJECT(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE,
                        (ObjectLabel::A, ObjectLabel::B, ObjectLabel::None))

#undef INSTANTIATE

#define FRAME(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data) \
  template struct ExcisionCenter<OBJECT(data), FRAME(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE,
                        (ObjectLabel::A, ObjectLabel::B, ObjectLabel::None),
                        (::Frame::Grid, ::Frame::Distorted))

#undef INSTANTIATE
#undef OBJECT
#undef FRAME
}  // namespace domain::Tags
