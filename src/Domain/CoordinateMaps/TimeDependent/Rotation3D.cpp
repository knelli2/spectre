// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/Rotation3D.hpp"

#include <cmath>
#include <ostream>
#include <pup.h>
#include <pup_stl.h>
#include <utility>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/TimeDependent/RotationMatrixHelpers.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Utilities/DereferenceWrapper.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdHelpers.hpp"

namespace domain::CoordinateMaps::TimeDependent {

Rotation3D::Rotation3D(std::string function_of_time_name) noexcept
    : f_of_t_name_(std::move(function_of_time_name)), stored_time_(0.0) {
  stored_matrix_and_deriv_[0] =
      Matrix{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  stored_matrix_and_deriv_[1] =
      Matrix{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
}

template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, 3> Rotation3D::operator()(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  if (not equal_within_roundoff(time, stored_time_)) {
    UpdateRotationMatrices<1>(make_not_null(&stored_matrix_and_deriv_), time,
                              functions_of_time, f_of_t_name_);
    stored_time_ = time;
  }

  std::array<tt::remove_cvref_wrap_t<T>, 3> result;
  for (size_t i = 0; i < 3; i++) {
    result[i] = stored_matrix_and_deriv_[0](i, 0) * source_coords[0];
    for (size_t j = 1; j < 3; j++) {
      result[i] += stored_matrix_and_deriv_[0](i, j) * source_coords[j];
    }
  }
  return result;
}

std::optional<std::array<double, 3>> Rotation3D::inverse(
    const std::array<double, 3>& target_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  if (not equal_within_roundoff(time, stored_time_)) {
    UpdateRotationMatrices<1>(make_not_null(&stored_matrix_and_deriv_), time,
                              functions_of_time, f_of_t_name_);
    stored_time_ = time;
  }

  // The inverse map uses the inverse rotation matrix, which is just the
  // transpose of the rotation matrix
  std::array<double, 3> result;
  for (size_t i = 0; i < 3; i++) {
    result[i] = stored_matrix_and_deriv_[0](0, i) * target_coords[0];
    for (size_t j = 1; j < 3; j++) {
      result[i] += stored_matrix_and_deriv_[0](j, i) * target_coords[j];
    }
  }
  return result;
}

template <typename T>
std::array<tt::remove_cvref_wrap_t<T>, 3> Rotation3D::frame_velocity(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  // The mapped coordinates (x,y,z) are related to the unmapped
  // coordinates (\xi, \eta, \zeta) by
  // (x,y,z) = dtR * (\xi, \eta, \zeta)
  // where dtR is the derivative of the rotation matrix
  if (not equal_within_roundoff(time, stored_time_)) {
    UpdateRotationMatrices<1>(make_not_null(&stored_matrix_and_deriv_), time,
                              functions_of_time, f_of_t_name_);
    stored_time_ = time;
  }

  std::array<tt::remove_cvref_wrap_t<T>, 3> result;
  for (size_t i = 0; i < 3; i++) {
    result[i] = stored_matrix_and_deriv_[1](i, 0) * source_coords[0];
    for (size_t j = 1; j < 3; j++) {
      result[i] += stored_matrix_and_deriv_[1](i, j) * source_coords[j];
    }
  }
  return result;
}

template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame> Rotation3D::jacobian(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  if (not equal_within_roundoff(time, stored_time_)) {
    UpdateRotationMatrices<1>(make_not_null(&stored_matrix_and_deriv_), time,
                              functions_of_time, f_of_t_name_);
    stored_time_ = time;
  }

  // Make tensor of zeros with correct type
  tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame> jacobian_matrix{
      make_with_value<tt::remove_cvref_wrap_t<T>>(
          dereference_wrapper(source_coords[0]), 0.0)};
  // The jacobian is just the rotation matrix
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      jacobian_matrix.get(i, j) = stored_matrix_and_deriv_[0](i, j);
    }
  }

  return jacobian_matrix;
}

template <typename T>
tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame>
Rotation3D::inv_jacobian(
    const std::array<T, 3>& source_coords, const double time,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) const noexcept {
  ASSERT(functions_of_time.find(f_of_t_name_) != functions_of_time.end(),
         "Could not find function of time: '"
             << f_of_t_name_ << "' in functions of time. Known functions are "
             << keys_of(functions_of_time));

  if (not equal_within_roundoff(time, stored_time_)) {
    UpdateRotationMatrices<1>(make_not_null(&stored_matrix_and_deriv_), time,
                              functions_of_time, f_of_t_name_);
    stored_time_ = time;
  }

  tnsr::Ij<tt::remove_cvref_wrap_t<T>, 3, Frame::NoFrame> inv_jacobian_matrix{
      make_with_value<tt::remove_cvref_wrap_t<T>>(
          dereference_wrapper(source_coords[0]), 0.0)};
  // The inverse jacobian is just the inverse rotation matrix, which is the
  // transpose of the rotation matrix.
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      inv_jacobian_matrix.get(i, j) = stored_matrix_and_deriv_[0](j, i);
    }
  }

  return inv_jacobian_matrix;
}

void Rotation3D::pup(PUP::er& p) noexcept {
  p | f_of_t_name_;
  p | stored_matrix_and_deriv_;
}

bool operator==(const Rotation3D& lhs, const Rotation3D& rhs) noexcept {
  return lhs.f_of_t_name_ == rhs.f_of_t_name_ and
         lhs.stored_matrix_and_deriv_ == rhs.stored_matrix_and_deriv_;
}

bool operator!=(const Rotation3D& lhs, const Rotation3D& rhs) noexcept {
  return not(lhs == rhs);
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE(_, data)                                                \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)>      \
  Rotation3D::operator()(                                                   \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const noexcept;                                \
  template std::array<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data)>      \
  Rotation3D::frame_velocity(                                               \
      const std::array<DTYPE(data), DIM(data)>& source_coords,              \
      const double time,                                                    \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const noexcept;                                \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),        \
                    Frame::NoFrame>                                         \
  Rotation3D::jacobian(                                                     \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const noexcept;                                \
  template tnsr::Ij<tt::remove_cvref_wrap_t<DTYPE(data)>, DIM(data),        \
                    Frame::NoFrame>                                         \
  Rotation3D::inv_jacobian(                                                 \
      const std::array<DTYPE(data), DIM(data)>& source_coords, double time, \
      const std::unordered_map<                                             \
          std::string,                                                      \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&        \
          functions_of_time) const noexcept;

GENERATE_INSTANTIATIONS(INSTANTIATE, (3),
                        (double, DataVector,
                         std::reference_wrapper<const double>,
                         std::reference_wrapper<const DataVector>))
#undef DIM
#undef DTYPE
#undef INSTANTIATE

}  // namespace domain::CoordinateMaps::TimeDependent
