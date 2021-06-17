// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <iostream>

template <typename T>
class Quaternion {
 public:
  Quaternion();
  explicit Quaternion(const T& a, const T& b, const T& c, const T& d);
  template <typename R>
  explicit Quaternion(const Quaternion<R>& quat_copy);
  template <typename R>
  explicit Quaternion(const std::array<R, 4>& array_copy);
  template <typename R>
  explicit Quaternion(const std::array<R, 3>& short_array_copy);

  Quaternion<T> real() const;
  Quaternion<T> unreal() const;
  Quaternion<T> conj() const;
  T norm() const;
  T normsqr() const;
  void normalize();

  template <typename R>
  T operator[](R i);
  template <typename R>
  T& operator[](R i);

  Quaternion<T>& operator=(const Quaternion<T>& rhs);
  template <typename R>
  Quaternion<T>& operator=(const Quaternion<R>& rhs);
  Quaternion<T>& operator=(const T& rhs);

  template <typename R, typename S>
  bool operator==(const Quaternion<R>& lhs, const Quternion<S>& rhs);
  template <typename R, typename S>
  bool operator==(const Quaternion<R>& lhs, const S& rhs);
  template <typename R, typename S>
  bool operator==(const R& lhs, const Quternion<S>& rhs);

  template <typename R, typename S>
  bool operator!=(const Quaternion<R>& lhs, const Quternion<S>& rhs);
  template <typename R, typename S>
  bool operator!=(const Quaternion<R>& lhs, const S& rhs);
  template <typename R, typename S>
  bool operator!=(const R& lhs, const Quternion<S>& rhs);

  template <typename R, typename S>
  Quaternion<T>& operator+(const Quaternion<R>& q1, const Quaternion<S>& q2);
  template <typename R, typename S>
  Quaternion<T>& operator+(const Quaternion<R>& q1, const S& rhs);
  template <typename R, typename S>
  Quaternion<T>& operator+(const R& lhs, const Quaternion<S>& q1);

  template <typename R>
  Quaternion<T>& operator+(const Quaternion<R>& q);

  template <typename R, typename S>
  Quaternion<T>& operator-(const Quaternion<R>& q1, const Quaternion<S>& q2);
  template <typename R, typename S>
  Quaternion<T>& operator-(const Quaternion<R>& q1, const S& rhs);
  template <typename R, typename S>
  Quaternion<T>& operator-(const R& lhs, const Quaternion<S>& q1);

  template <typename R>
  Quaternion<T>& operator-(const Quaternion<R>& q);

  template <typename R, typename S>
  Quaternion<T>& operator*(const Quaternion<R>& q1, const Quaternion<S>& q2);
  template <typename R, typename S>
  Quaternion<T>& operator*(const Quaternion<R>& q1, const S& rhs);
  template <typename R, typename S>
  Quaternion<T>& operator*(const R& lhs, const Quaternion<S>& q1);

  template <typename R, typename S>
  Quaternion<T>& operator/(const Quaternion<R>& q1, const Quaternion<S>& q2);
  template <typename R, typename S>
  Quaternion<T>& operator/(const Quaternion<R>& q1, const S& rhs);
  template <typename R, typename S>
  Quaternion<T>& operator/(const R& lhs, const Quaternion<S>& q1);

  template <typename R>
  Quaternion<T>& operator+=(const Quaternion<R>& rhs);
  template <typename R>
  Quaternion<T>& operator+=(R& rhs);

  template <typename R>
  Quaternion<T>& operator-=(const Quaternion<R>& rhs);
  template <typename R>
  Quaternion<T>& operator-=(R& rhs);

  template <typename R>
  Quaternion<T>& operator*=(const Quaternion<R>& rhs);
  template <typename R>
  Quaternion<T>& operator*=(R& rhs);

  template <typename R>
  Quaternion<T>& operator/=(const Quaternion<R>& rhs);
  template <typename R>
  Quaternion<T>& operator/=(R& rhs);

  ostream& operator << (ostream& os, const Quaternion<T> q);

 private:
  a_;
  b_;
  c_;
  d_;
};
