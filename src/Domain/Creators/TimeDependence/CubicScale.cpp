// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependence/CubicScale.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/MapInstantiationMacros.hpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"

namespace domain {
namespace creators::time_dependence {
template <size_t MeshDim>
CubicScale<MeshDim>::CubicScale(const double initial_time,
                                const double outer_boundary,
                                bool use_linear_scaling,
                                const std::array<double, 2>& initial_expansion,
                                const std::array<double, 2>& velocity,
                                const std::array<double, 2>& acceleration)
    : initial_time_(initial_time),
      outer_boundary_(outer_boundary),
      use_linear_scaling_(use_linear_scaling),
      initial_expansion_(initial_expansion),
      velocity_(velocity),
      acceleration_(acceleration) {
  // This makes the function name unique because this function of time doesn't
  // expire. This also encodes initial data info in the name for diagnostic
  // purposes
  if (use_linear_scaling_) {
    // If we are using linear scaling, then the names must be the same so it
    // doesn't matter about initial conditions
    functions_of_time_names_[0] = "CubicScale";
    functions_of_time_names_[1] = "CubicScale";
  } else {
    functions_of_time_names_[0] =
        "CubicScaleA::a=" + get_output(initial_expansion_[0]) +
        "::dta=" + get_output(velocity_[0]) +
        "::d2ta=" + get_output(acceleration_[0]) +
        "::t_0=" + get_output(initial_time_) +
        "::outer_boundary=" + get_output(outer_boundary_);
    functions_of_time_names_[1] =
        "CubicScaleB::a=" + get_output(initial_expansion_[1]) +
        "::dta=" + get_output(velocity_[1]) +
        "::d2ta=" + get_output(acceleration_[1]) +
        "::t_0=" + get_output(initial_time_) +
        "::outer_boundary=" + get_output(outer_boundary_);
  }
}

template <size_t MeshDim>
std::unique_ptr<TimeDependence<MeshDim>> CubicScale<MeshDim>::get_clone()
    const {
  return std::make_unique<CubicScale>(initial_time_, outer_boundary_,
                                      use_linear_scaling_, initial_expansion_,
                                      velocity_, acceleration_);
}

template <size_t MeshDim>
std::vector<std::unique_ptr<
    domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, MeshDim>>>
CubicScale<MeshDim>::block_maps(const size_t number_of_blocks) const {
  ASSERT(number_of_blocks > 0, "Must have at least one block to create.");
  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, MeshDim>>>
      result{number_of_blocks};
  result[0] = std::make_unique<MapForComposition>(map_for_composition());
  for (size_t i = 1; i < number_of_blocks; ++i) {
    result[i] = result[0]->get_clone();
  }
  return result;
}

template <size_t MeshDim>
std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
CubicScale<MeshDim>::functions_of_time(
    const std::vector<std::pair<std::string, double>>& initial_expiration_times)
    const {
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      result{};

  // Use a 3rd deriv function of time so that it can be used with a control
  // system.
  if (initial_expiration_times.size()) {
    // initial_expiration_times is ordered such that object A is first, and
    // object B is second
    ASSERT(initial_expiration_times.size() == 2,
           "There are only 2 functions of time for a CubicScale time "
           "dependence, however, " +
               get_output(initial_expiration_times.size()) +
               " initial exiration times were supplied.");
    for (size_t i = 0; i < initial_expiration_times.size(); i++) {
      functions_of_time_names_[i] = initial_expiration_times[i].first;
      expiration_times_[i] = initial_expiration_times[i].second;
    }
  }

  // If we are using linear scaling, the function of time names will be the same
  // so this first assignment will be overwritten by the second assignmnet.
  // This is expected.
  result[functions_of_time_names_[0]] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<3>>(
          initial_time_,
          std::array<DataVector, 4>{{{initial_expansion_[0]},
                                     {velocity_[0]},
                                     {acceleration_[0]},
                                     {0.0}}},
          expiration_times_[0]);
  result[functions_of_time_names_[1]] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<3>>(
          initial_time_,
          std::array<DataVector, 4>{{{initial_expansion_[1]},
                                     {velocity_[1]},
                                     {acceleration_[1]},
                                     {0.0}}},
          expiration_times_[0]);
  return result;
}

template <size_t MeshDim>
auto CubicScale<MeshDim>::map_for_composition() const -> MapForComposition {
  return MapForComposition{CubicScaleMap{outer_boundary_,
                                         functions_of_time_names_[0],
                                         functions_of_time_names_[1]}};
}

template <size_t Dim>
bool operator==(const CubicScale<Dim>& lhs, const CubicScale<Dim>& rhs) {
  return lhs.initial_time_ == rhs.initial_time_ and
         lhs.outer_boundary_ == rhs.outer_boundary_ and
         lhs.functions_of_time_names_ == rhs.functions_of_time_names_ and
         lhs.initial_expansion_ == rhs.initial_expansion_ and
         lhs.velocity_ == rhs.velocity_ and
         lhs.acceleration_ == rhs.acceleration_;
}

template <size_t Dim>
bool operator!=(const CubicScale<Dim>& lhs, const CubicScale<Dim>& rhs) {
  return not(lhs == rhs);
}

#define GET_DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                               \
  template class CubicScale<GET_DIM(data)>;                                  \
  template bool operator==<GET_DIM(data)>(const CubicScale<GET_DIM(data)>&,  \
                                          const CubicScale<GET_DIM(data)>&); \
  template bool operator!=<GET_DIM(data)>(const CubicScale<GET_DIM(data)>&,  \
                                          const CubicScale<GET_DIM(data)>&);

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef GET_DIM
}  // namespace creators::time_dependence

INSTANTIATE_MAPS_FUNCTIONS(((CoordinateMaps::TimeDependent::CubicScale<1>),
                            (CoordinateMaps::TimeDependent::CubicScale<2>),
                            (CoordinateMaps::TimeDependent::CubicScale<3>)),
                           (Frame::Grid), (Frame::Inertial),
                           (double, DataVector))
}  // namespace domain
