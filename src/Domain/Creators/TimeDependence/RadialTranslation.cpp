// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependence/RadialTranslation.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/MapInstantiationMacros.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace domain {
namespace creators::time_dependence {

template <size_t MeshDim>
RadialTranslation<MeshDim>::RadialTranslation(
    double initial_time, std::vector<double> inner_map_initial_values,
    const std::optional<OuterMapOptions>& outer_map_options,
    const Options::Context& context)
    : initial_time_(initial_time),
      inner_initial_values_(std::move(inner_map_initial_values)),
      outer_map_options_(outer_map_options) {
  if (outer_map_options_.has_value()) {
    if (outer_map_options_->inner_outer_radius[0] >
        outer_map_options_->inner_outer_radius[1]) {
      PARSE_ERROR(context, "Inner radius "
                               << outer_map_options_->inner_outer_radius[0]
                               << " must be smaller than the outer radius "
                               << outer_map_options_->inner_outer_radius[1]);
    }
  }
}

template <size_t MeshDim>
std::unique_ptr<TimeDependence<MeshDim>> RadialTranslation<MeshDim>::get_clone()
    const {
  return std::make_unique<RadialTranslation>(
      initial_time_, inner_initial_values_, outer_map_options_);
}

template <size_t MeshDim>
std::vector<std::unique_ptr<
    domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, MeshDim>>>
RadialTranslation<MeshDim>::block_maps_grid_to_inertial(
    const size_t number_of_blocks) const {
  ASSERT(number_of_blocks > 0, "Must have at least one block to create.");
  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, MeshDim>>>
      result{number_of_blocks};
  result[0] = std::make_unique<GridToInertialMap>(grid_to_inertial_map());
  for (size_t i = 1; i < number_of_blocks; ++i) {
    result[i] = result[0]->get_clone();
  }
  return result;
}

template <size_t MeshDim>
std::vector<std::unique_ptr<
    domain::CoordinateMapBase<Frame::Grid, Frame::Distorted, MeshDim>>>
RadialTranslation<MeshDim>::block_maps_grid_to_distorted(
    const size_t number_of_blocks) const {
  ASSERT(number_of_blocks > 0, "Must have at least one block to create.");
  using ptr_type =
      domain::CoordinateMapBase<Frame::Grid, Frame::Distorted, MeshDim>;
  std::vector<std::unique_ptr<ptr_type>> result{number_of_blocks};

  return result;
}

template <size_t MeshDim>
std::vector<std::unique_ptr<
    domain::CoordinateMapBase<Frame::Distorted, Frame::Inertial, MeshDim>>>
RadialTranslation<MeshDim>::block_maps_distorted_to_inertial(
    const size_t number_of_blocks) const {
  ASSERT(number_of_blocks > 0, "Must have at least one block to create.");
  using ptr_type =
      domain::CoordinateMapBase<Frame::Distorted, Frame::Inertial, MeshDim>;
  std::vector<std::unique_ptr<ptr_type>> result{number_of_blocks};

  return result;
}

template <size_t MeshDim>
std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
RadialTranslation<MeshDim>::functions_of_time(
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const {
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      result{};

  // We use a `PiecewisePolynomial` with 2 derivs since some transformations
  // between different frames for moving meshes can require Hessians.

  // Functions of time don't expire by default
  double expiration_time = std::numeric_limits<double>::infinity();

  // If we have control systems, overwrite the expiration time with the one
  // supplied by the control system
  if (initial_expiration_times.count(function_of_time_name_) == 1) {
    expiration_time = initial_expiration_times.at(function_of_time_name_);
  }

  std::array<DataVector, 3> initial_values = make_array<3>(
      DataVector{outer_map_options_.has_value() ? 2_st : 1_st, 0.0});

  for (size_t i = 0; i < inner_initial_values_.size(); i++) {
    gsl::at(initial_values, i)[0] = inner_initial_values_[i];
  }

  if (outer_map_options_.has_value()) {
    for (size_t i = 0; i < outer_map_options_->initial_values.size(); i++) {
      gsl::at(initial_values, i)[1] = outer_map_options_->initial_values[i];
    }
  }

  result[function_of_time_name_] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
          initial_time_, std::move(initial_values), expiration_time);

  return result;
}

template <size_t MeshDim>
auto RadialTranslation<MeshDim>::grid_to_inertial_map() const
    -> GridToInertialMap {
  if (outer_map_options_.has_value()) {
    return GridToInertialMap{
        RadialTranslationMap{function_of_time_name_,
                             make_array<MeshDim>(0.0),
                             {outer_map_options_->inner_outer_radius}}};
  } else {
    return GridToInertialMap{
        RadialTranslationMap{function_of_time_name_, make_array<MeshDim>(0.0)}};
  }
}

template <size_t Dim>
bool operator==(const RadialTranslation<Dim>& lhs,
                const RadialTranslation<Dim>& rhs) {
  // If distorted_and_inertial_frames_are_equal_, we short-circuit
  // evaluating velocity_distorted_to_inertial_ below.
  return lhs.initial_time_ == rhs.initial_time_ and
         lhs.inner_initial_values_ == rhs.inner_initial_values_ and
         lhs.outer_map_options_.has_value() ==
             rhs.outer_map_options_.has_value() and
         (lhs.outer_map_options_.has_value() and
          lhs.outer_map_options_->initial_values ==
              rhs.outer_map_options_->initial_values and
          lhs.outer_map_options_->inner_outer_radius ==
              rhs.outer_map_options_->inner_outer_radius);
}

template <size_t Dim>
bool operator!=(const RadialTranslation<Dim>& lhs,
                const RadialTranslation<Dim>& rhs) {
  return not(lhs == rhs);
}

#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                  \
  template class RadialTranslation<GET_DIM(data)>;              \
  template bool operator==                                      \
      <GET_DIM(data)>(const RadialTranslation<GET_DIM(data)>&,  \
                      const RadialTranslation<GET_DIM(data)>&); \
  template bool operator!=                                      \
      <GET_DIM(data)>(const RadialTranslation<GET_DIM(data)>&,  \
                      const RadialTranslation<GET_DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATION, (2, 3))

#undef GET_DIM
#undef INSTANTIATION
}  // namespace creators::time_dependence

template <size_t MeshDim>
using RadialTranslation =
    CoordinateMaps::TimeDependent::RadialTranslation<MeshDim>;

INSTANTIATE_MAPS_FUNCTIONS(((RadialTranslation<2>), (RadialTranslation<3>)),
                           (Frame::Grid), (Frame::Distorted, Frame::Inertial),
                           (double, DataVector))
INSTANTIATE_MAPS_FUNCTIONS(((RadialTranslation<2>), (RadialTranslation<3>)),
                           (Frame::Distorted), (Frame::Inertial),
                           (double, DataVector))

}  // namespace domain
