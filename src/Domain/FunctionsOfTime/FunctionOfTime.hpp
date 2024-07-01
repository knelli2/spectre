// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <memory>
#include <pup.h>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"

namespace domain {
/// \ingroup ComputationalDomainGroup
/// \brief Contains functions of time to support the dual frame system.
namespace FunctionsOfTime {
/// \ingroup ComputationalDomainGroup
/// \brief Base class for FunctionsOfTime
///
/// A FunctionOfTime is a function that will return the same value for
/// a time `t`, regardless of when that function is called during a run
/// (provided that the time `t` is in the domain of validity of the function).
/// All FunctionsOfTime have members
///   - `func`, that returns a `std::array<DataVector, 1>`
///   - `func_and_deriv`, that returns a `std::array<DataVector, 2>`
///   - `func_and_2_derivs`, that returns a `std::array<DataVector, 3>`
///
/// The DataVectors that are returned can be of any size: e.g. a scalar
/// FunctionOfTime will have DataVectors with one component and a 3-vector
/// FunctionOfTime will have DataVectors with three components.
///
/// The domain of validity of the function is given by the `time_bounds` member
/// function.
///
/// The function and all of its derivatives are left-continuous, that
/// is, when evaluated at a time when the function was updated, they
/// return the values just before the update, ignoring the updated
/// value.
class FunctionOfTime : public PUP::able {
 public:
  FunctionOfTime() = default;
  FunctionOfTime(FunctionOfTime&&) = default;
  FunctionOfTime& operator=(FunctionOfTime&&) = default;
  FunctionOfTime(const FunctionOfTime&) = default;
  FunctionOfTime& operator=(const FunctionOfTime&) = default;
  ~FunctionOfTime() override = default;

  virtual auto get_clone() const -> std::unique_ptr<FunctionOfTime> = 0;

  /// Returns the domain of validity of the function.
  /// For FunctionsOfTime that allow a small amount of time extrapolation,
  /// `time_bounds` tells you the bounds including the allowed extrapolation
  /// interval.
  virtual std::array<double, 2> time_bounds() const = 0;

  /// \brief The first expiration time after \p time.
  ///
  /// \details For non-updatable functions, this returns infinity. For
  /// updatable functions, the first expiration time after \p time is
  /// found by determining the update immediately before \p time. The
  /// expiration time of this update is what is returned. If \p time
  /// happens to be an update itself, then the expiration of that
  /// update is returned.
  virtual double expiration_after(double time) const = 0;

  /// Updates the maximum derivative of the FunctionOfTime at a given time while
  /// also resetting the expiration. By default, a FunctionOfTime cannot be
  /// updated.
  virtual void update(double /*time_of_update*/,
                      DataVector /*updated_max_deriv*/,
                      double /*next_expiration_time*/) {
    ERROR("Cannot update this FunctionOfTime.");
  }

  /// The DataVector can be of any size
  virtual std::array<DataVector, 1> func(double t) const = 0;
  /// The DataVector can be of any size
  virtual std::array<DataVector, 2> func_and_deriv(double t) const = 0;
  /// The DataVector can be of any size
  virtual std::array<DataVector, 3> func_and_2_derivs(double t) const = 0;

  virtual std::vector<DataVector> func_and_all_derivs(double t) const {
    std::array<DataVector, 3> tmp_func_and_2_derivs = func_and_2_derivs(t);

    return std::vector{std::move(tmp_func_and_2_derivs[0]),
                       std::move(tmp_func_and_2_derivs[1]),
                       std::move(tmp_func_and_2_derivs[2])};
  }

  WRAPPED_PUPable_abstract(FunctionOfTime);  // NOLINT
};
}  // namespace FunctionsOfTime

using FunctionsOfTimeMap = std::unordered_map<
    std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>;

}  // namespace domain
