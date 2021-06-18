// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>

#include "Utilities/Gsl.hpp"

template <typename T>
class Quaternion {
 public:
  Quaternion() noexcept
      : data_(std::array<T, 4>{static_cast<T>(0), static_cast<T>(0),
                               static_cast<T>(0), static_cast<T>(0)}) {}
  Quaternion(const T& a, const T& b, const T& c, const T& d) noexcept
      : data_(std::array<T, 4>{a, b, c, d}) {}
  Quaternion(const Quaternion<T>& quat_copy) noexcept
      : data_(quat_copy.array()) {}
  template <typename R>
  Quaternion(const Quaternion<R>& quat_copy) noexcept
      : data_(std::array<T, 4>{
            static_cast<T>(quat_copy[0]), static_cast<T>(quat_copy[1]),
            static_cast<T>(quat_copy[2]), static_cast<T>(quat_copy[3])}) {}
  template <typename R>
  Quaternion(const std::array<R, 4>& array_copy) noexcept
      : data_(std::array<T, 4>{
            static_cast<T>(array_copy[0]), static_cast<T>(array_copy[1]),
            static_cast<T>(array_copy[2]), static_cast<T>(array_copy[3])}) {}
  template <typename R>
  Quaternion(const std::array<R, 3>& short_array_copy) noexcept
      : data_(std::array<T, 4>{static_cast<T>(0),
                               static_cast<T>(short_array_copy[0]),
                               static_cast<T>(short_array_copy[1]),
                               static_cast<T>(short_array_copy[2])}) {}

  Quaternion<T> real() const noexcept {
    return Quaternion<T>(data_[0], static_cast<T>(0), static_cast<T>(0),
                         static_cast<T>(0));
  }
  Quaternion<T> unreal() const noexcept {
    return Quaternion<T>(static_cast<T>(0), data_[1], data_[2], data_[3]);
  }
  Quaternion<T> conj() const noexcept {
    return Quaternion<T>(data_[0], -data_[1], -data_[2], -data_[3]);
  }
  T norm() const noexcept {
    return static_cast<T>(sqrt(data_[0] * data_[0] + data_[1] * data_[1] +
                               data_[2] * data_[2] + data_[3] * data_[3]));
  }
  T normsqr() const noexcept {
    return static_cast<T>(data_[0] * data_[0] + data_[1] * data_[1] +
                          data_[2] * data_[2] + data_[3] * data_[3]);
  }
  Quaternion<T> inv() const noexcept { return conj() / normsqr(); }
  void unit() noexcept {
    T norm = norm();
    for (auto elem : data_) {
      elem /= norm;
    }
  }

  using iterator = typename std::array<T, 4>::iterator;
  using const_iterator = typename std::array<T, 4>::const_iterator;

  iterator begin() noexcept { return data_.begin(); }
  const_iterator begin() const noexcept { return data_.begin(); }

  iterator end() noexcept { return data_.end(); }
  const_iterator end() const noexcept { return data_.end(); }

  std::array<T, 4> array() const noexcept { return data_; }

  T& operator[](const size_t i) noexcept { return gsl::at(data_, i); }
  const T& operator[](const size_t i) const noexcept {
    return gsl::at(data_, i);
  }

  Quaternion<T>& operator=(const Quaternion<T>& rhs) {
    data_ = rhs.array();
    return *this;
  }
  template <typename R>
  Quaternion<T>& operator=(const Quaternion<R>& rhs) {
    data_[0] = static_cast<T>(rhs[0]);
    data_[1] = static_cast<T>(rhs[1]);
    data_[2] = static_cast<T>(rhs[2]);
    data_[3] = static_cast<T>(rhs[3]);
    return *this;
  }
  Quaternion<T>& operator=(const T& rhs) {
    data_[0] = rhs;
    data_[1] = static_cast<T>(0);
    data_[2] = static_cast<T>(0);
    data_[3] = static_cast<T>(0);
    return *this;
  }
  template <typename R>
  Quaternion<T>& operator=(const R& rhs) {
    data_[0] = static_cast<T>(rhs);
    data_[1] = static_cast<T>(0);
    data_[2] = static_cast<T>(0);
    data_[3] = static_cast<T>(0);
    return *this;
  }

  bool operator==(const Quaternion<T>& rhs) const {
    return data_ == rhs.array();
  }
  bool operator==(const T& rhs) const {
    return data_ == std::array<T, 4>{rhs, static_cast<T>(0), static_cast<T>(0),
                                     static_cast<T>(0)};
  }

  bool operator!=(const Quaternion<T>& rhs) const { return not(*this == rhs); }
  bool operator!=(const T& rhs) const { return not(*this == rhs); }

  Quaternion<T>& operator+=(const Quaternion<T>& rhs) {
    auto rhs_it = rhs.begin();
    auto data_it = data_.begin();
    for (; data_it != data_.end(); rhs_it++, data_it++) {
      *data_it += *rhs_it;
    }
    return *this;
  }
  Quaternion<T>& operator+=(const T& rhs) {
    data_[0] += rhs;
    return *this;
  }

  Quaternion<T>& operator-=(const Quaternion<T>& rhs) {
    auto rhs_it = rhs.begin();
    auto data_it = data_.begin();
    for (; data_it != data_.end(); rhs_it++, data_it++) {
      *data_it -= *rhs_it;
    }
    return *this;
  }
  Quaternion<T>& operator-=(T& rhs) {
    data_[0] -= rhs;
    return *this;
  }

  Quaternion<T>& operator*=(const Quaternion<T>& rhs) {
    Quaternion<T> tmpq = *this * rhs;
    data_[0] = tmpq[0];
    data_[1] = tmpq[1];
    data_[2] = tmpq[2];
    data_[3] = tmpq[3];
    return *this;
  }
  Quaternion<T>& operator*=(const T& rhs) {
    for (auto elem : data_) {
      elem *= rhs;
    }
    return *this;
  }

  Quaternion<T>& operator/=(const Quaternion<T>& rhs) {
    Quaternion<T> tmpq = *this * rhs.inv();
    data_[0] = tmpq[0];
    data_[1] = tmpq[1];
    data_[2] = tmpq[2];
    data_[3] = tmpq[3];
    return *this;
  }
  Quaternion<T>& operator/=(T& rhs) {
    for (auto elem : data_) {
      elem /= rhs;
    }
    return *this;
  }

 private:
  std::array<T, 4> data_;
};

template <typename T>
Quaternion<T> operator+(const Quaternion<T>& q1, const Quaternion<T>& q2) {
  return Quaternion<T>(q1[0] + q2[0], q1[1] + q2[1], q1[2] + q2[2],
                       q1[3] + q2[3]);
}
template <typename T>
Quaternion<T> operator+(const Quaternion<T>& q, const T& rhs) {
  return Quaternion<T>(q[0] + rhs, q[1], q[2], q[3]);
}
template <typename T>
Quaternion<T> operator+(const T& lhs, const Quaternion<T>& q) {
  return Quaternion<T>(lhs + q[0], q[1], q[2], q[3]);
}

template <typename T>
Quaternion<T>& operator+(const Quaternion<T>& q) {
  return q;
}

template <typename T>
Quaternion<T> operator-(const Quaternion<T>& q1, const Quaternion<T>& q2) {
  return Quaternion<T>(q1[0] - q2[0], q1[1] - q2[1], q1[2] - q2[2],
                       q1[3] - q2[3]);
}
template <typename T>
Quaternion<T> operator-(const Quaternion<T>& q, const T& rhs) {
  return Quaternion<T>(q[0] - rhs, q[1], q[2], q[3]);
}
template <typename T>
Quaternion<T> operator-(const T& lhs, const Quaternion<T>& q) {
  return Quaternion<T>(lhs - q[0], -q[1], -q[2], -q[3]);
}

template <typename T>
Quaternion<T> operator-(const Quaternion<T>& q) {
  return Quaternion<T>(-q[0], -q[1], -q[2], -q[3]);
}

template <typename T>
Quaternion<T> operator*(const Quaternion<T>& q1, const Quaternion<T>& q2) {
  return Quaternion<T>(
      q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3],
      q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2],
      q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1],
      q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]);
}
template <typename T>
Quaternion<T> operator*(const Quaternion<T>& q, const T& rhs) {
  return Quaternion<T>(q[0] * rhs, q[1] * rhs, q[2] * rhs, q[3] * rhs);
}
template <typename T>
Quaternion<T> operator*(const T& lhs, const Quaternion<T>& q) {
  return Quaternion<T>(q[0] * lhs, q[1] * lhs, q[2] * lhs, q[3] * lhs);
}

template <typename T>
Quaternion<T> operator/(const Quaternion<T>& q1, const Quaternion<T>& q2) {
  return q1 * q2.inv();
}
template <typename T>
Quaternion<T> operator/(const Quaternion<T>& q, const T& rhs) {
  return Quaternion<T>(q[0] / rhs, q[1] / rhs, q[2] / rhs, q[3] / rhs);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Quaternion<T> q) {
  os << "(" << q[0] << "," << q[1] << "," << q[2] << "," << q[3] << ")";
  return os;
}
