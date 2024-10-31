// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/DiscontinuousGalerkin/ElementState.hpp"

#include <string>

#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

std::ostream& operator<<(std::ostream& os, const ElementState ordering) {
  switch (ordering) {
    case ElementState::ChuggingAlong:
      return os << "ChuggingAlong";
    case ElementState::WaitingForBoundaryDataInAction:
      return os << "WaitingForBoundaryDataInAction";
    case ElementState::WaitingForBoundaryDataInPostProcessor:
      return os << "WaitingForBoundaryDataInPostProcessor";
    case ElementState::WaitingForFoTInAction:
      return os << "WaitingForFoTInAction";
    case ElementState::WaitingForFoTInPostProcessor:
      return os << "WaitingForFoTInPostProcessor";
    case ElementState::WaitingForMeasurementTimescales1:
      return os << "WaitingForMeasurementTimescales";
    default:
      ERROR("Unknown ElementState");
  }
}

template <>
ElementState Options::create_from_yaml<ElementState>::create<void>(
    const Options::Option& options) {
  const auto ordering = options.parse_as<std::string>();
  if (ordering == "ChuggingAlong") {
    return ElementState::ChuggingAlong;
  } else if (ordering == "WaitingForBoundaryDataInAction") {
    return ElementState::WaitingForBoundaryDataInAction;
  } else if (ordering == "WaitingForBoundaryDataInPostProcessor") {
    return ElementState::WaitingForBoundaryDataInPostProcessor;
  } else if (ordering == "WaitingForFoTInAction") {
    return ElementState::WaitingForFoTInAction;
  } else if (ordering == "WaitingForFoTInPostProcessor") {
    return ElementState::WaitingForFoTInPostProcessor;
  } else if (ordering == "WaitingForMeasurementTimescales") {
    return ElementState::WaitingForMeasurementTimescales1;
  }
  PARSE_ERROR(options.context(), "ElementState wrong");
}
