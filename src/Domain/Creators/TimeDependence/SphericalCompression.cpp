// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependence/SphericalCompression.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/MapInstantiationMacros.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Options/Options.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace domain {
namespace creators::time_dependence {

SphericalCompression::SphericalCompression(
    const double initial_time, const double min_radius, const double max_radius,
    const std::array<double, 3> center, const double initial_value,
    const double initial_velocity, const double initial_acceleration,
    const Options::Context& context)
    : initial_time_(initial_time),
      min_radius_(min_radius),
      max_radius_(max_radius),
      center_(center),
      initial_value_(initial_value),
      initial_velocity_(initial_velocity),
      initial_acceleration_(initial_acceleration) {
  if (min_radius >= max_radius) {
    PARSE_ERROR(context,
                "Tried to create a SphericalCompression TimeDependence, but "
                "the minimum radius ("
                    << min_radius << ") is not less than the maximum radius ("
                    << max_radius << ")");
  }
  // This makes the function name unique because this function of time doesn't
  // expire. This also encodes all initial data info in the name for diagnostic
  // purposes
  // clang-format off
  function_of_time_name_ =
      "SpherialCompression"s +
      "::r_min="s + get_output(min_radius_) +
      "::r_max="s + get_output(max_radius_) +
      "::center="s + get_output(center_) +
      "::value="s + get_output(initial_value_) +
      "::dtvalue="s + get_output(initial_velocity_) +
      "::d2tvalue="s + get_output(initial_acceleration_) +
      "::t_0="s + get_output(initial_time_);
  // clang-format on
}

std::unique_ptr<TimeDependence<3>> SphericalCompression::get_clone() const {
  return std::make_unique<SphericalCompression>(
      initial_time_, min_radius_, max_radius_, center_, initial_value_,
      initial_velocity_, initial_acceleration_);
}

std::vector<
    std::unique_ptr<domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>>>
SphericalCompression::block_maps(const size_t number_of_blocks) const {
  ASSERT(number_of_blocks > 0,
         "Must have at least one block on which to create a map.");
  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>>>
      result{number_of_blocks};
  result[0] = std::make_unique<MapForComposition>(map_for_composition());
  for (size_t i = 1; i < number_of_blocks; ++i) {
    result[i] = result[0]->get_clone();
  }
  return result;
}

std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
SphericalCompression::functions_of_time(
    const std::vector<std::pair<std::string, double>>& initial_expiration_times)
    const {
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      result{};

  if (initial_expiration_times.size()) {
    ASSERT(initial_expiration_times.size() == 1,
           "There is only 1 function of time for a SphericalCompression time "
           "dependence, however, " +
               get_output(initial_expiration_times.size()) +
               " initial exiration times were supplied.");
    // initial_expiration_times should only have one item so this is ok
    function_of_time_name_ = initial_expiration_times[0].first;
    expiration_time_ = initial_expiration_times[0].second;
  }

  result[function_of_time_name_] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<3>>(
          initial_time_,
          std::array<DataVector, 4>{{{initial_value_},
                                     {initial_velocity_},
                                     {initial_acceleration_},
                                     {0.0}}},
          expiration_time_);
  return result;
}

auto SphericalCompression::map_for_composition() const -> MapForComposition {
  return MapForComposition{SphericalCompressionMap{
      function_of_time_name_, min_radius_, max_radius_, center_}};
}

bool operator==(const SphericalCompression& lhs,
                const SphericalCompression& rhs) {
  return lhs.initial_time_ == rhs.initial_time_ and
         lhs.min_radius_ == rhs.min_radius_ and
         lhs.max_radius_ == rhs.max_radius_ and lhs.center_ == rhs.center_ and
         lhs.initial_value_ == rhs.initial_value_ and
         lhs.initial_velocity_ == rhs.initial_velocity_ and
         lhs.initial_acceleration_ == rhs.initial_acceleration_ and
         lhs.function_of_time_name_ == rhs.function_of_time_name_;
}

bool operator!=(const SphericalCompression& lhs,
                const SphericalCompression& rhs) {
  return not(lhs == rhs);
}
}  // namespace creators::time_dependence

using SphericalCompressionMap3d =
    CoordinateMaps::TimeDependent::SphericalCompression<false>;

INSTANTIATE_MAPS_FUNCTIONS(((SphericalCompressionMap3d)), (Frame::Grid),
                           (Frame::Inertial), (double, DataVector))

}  // namespace domain
