// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>

#include "DataStructures/Quaternion.hpp"

SPECTRE_TEST_CASE("Unit.DataStructures.Quaternion", "[Unit]") {
  {
    INFO("Check constructors and = operator (and also implicitly .array() )");
    Quaternion<double> quat1;
    Quaternion<double> quat2(1.0, 2.0, 3.0, 4.0);
    Quaternion<double> quat3(quat2);
    Quaternion<double> quat4(std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    Quaternion<double> quat5(std::array<double, 3>({1.0, 2.0, 3.0}));
    Quaternion<int> quat6(1, 2, 3, 4);
    Quaternion<double> quat7(quat6);
    Quaternion<double> quat8(std::array<int, 4>({1, 2, 3, 4}));
    Quaternion<double> quat9(std::array<int, 3>({1, 2, 3}));
    Quaternion<double> quat10(8.8);
    Quaternion<double> quat11(7);
    CHECK(quat1.array() == std::array<double, 4>({0.0, 0.0, 0.0, 0.0}));
    CHECK(quat2.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat3.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat4.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat5.array() == std::array<double, 4>({0.0, 1.0, 2.0, 3.0}));
    CHECK(quat6.array() == std::array<int, 4>({1, 2, 3, 4}));
    CHECK(quat7.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat8.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat9.array() == std::array<double, 4>({0.0, 1.0, 2.0, 3.0}));
    CHECK(quat10.array() == std::array<double, 4>({8.8, 0.0, 0.0, 0.0}));
    CHECK(quat11.array() == std::array<double, 4>({7.0, 0.0, 0.0, 0.0}));

    quat1 = quat2;  // assignment of quaternion with same template type
    quat3 = quat6;  // assignment of quaternion with different template type
    quat4 = 8.8;    // assignment of literal with same template type
    quat5 = 7;      // assignment of literal with different template type
    // assignment of std::array with same template type
    quat7 = std::array<double, 4>({1.0, 2.0, 3.0, 4.0});
    // assignment of std::array with different template type
    quat8 = std::array<int, 4>({1, 2, 3, 4});
    CHECK(quat1.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat3.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat4.array() == std::array<double, 4>({8.8, 0.0, 0.0, 0.0}));
    CHECK(quat5.array() == std::array<double, 4>({7.0, 0.0, 0.0, 0.0}));
    CHECK(quat7.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
    CHECK(quat8.array() == std::array<double, 4>({1.0, 2.0, 3.0, 4.0}));
  }

  {
    INFO("Check member functions and ==, != operators");
    Approx custom_approx = Approx::custom().epsilon(1.e-15);
    Quaternion<double> quat(1.0, 2.0, 3.0, 4.0);
    CHECK(quat.real() == Quaternion<double>(1.0, 0.0, 0.0, 0.0));
    CHECK(quat.unreal() == Quaternion<double>(0.0, 2.0, 3.0, 4.0));
    CHECK(quat.conj() == Quaternion<double>(1.0, -2.0, -3.0, -4.0));
    CHECK(quat.norm() == custom_approx(sqrt(30.0)));
    CHECK(quat.normsqr() == custom_approx(30.0));
    CHECK(quat.inv() == Quaternion<double>(1.0 / 30.0, -2.0 / 30.0, -3.0 / 30.0,
                                           -4.0 / 30.0));
    quat.unit();
    CHECK(quat.array() ==
          std::array<double, 4>({1.0 / sqrt(30.0), 2.0 / sqrt(30.0),
                                 3.0 / sqrt(30.0), 4.0 / sqrt(30.0)}));
    CHECK(quat.norm() == custom_approx(1.0));

    Quaternion<double> quat2(4.3, 0.0, 0.0, 0.0);
    // cannot use custom_approx here because compiler complains about casting
    // between Catch::Detail::Approx to double
    CHECK(quat2 == 4.3);
    CHECK_FALSE(quat == quat2);
    CHECK(quat != quat2);
    CHECK_FALSE(quat == 1.0);
    CHECK(quat != 1.0);
  }

  {
    INFO("Check iterators and [] operator")
    Quaternion<double> quat(1.1, 1.2, 1.3, 1.4);
    auto it = quat.begin();
    auto it_end = quat.end();
    size_t i = 0;
    for (; it != it_end; it++, i++) {
      CHECK(*it == approx(quat[i]));
    }
  }

  {
    INFO("Check + - * / operators")
    Quaternion<double> quat1(1.1, 1.2, 1.3, 1.4);
    Quaternion<double> quat2(-9.8, -9.7, 9.6, -9.5);
    Quaternion<double> quat3;
    Quaternion<double> quat4;
    Quaternion<double> quat5;
    Quaternion<double> quat6;

    quat3 = quat1 + quat2;
    quat4 = quat2 + quat1;
    quat5 = quat1 + 2.5;
    quat6 = 2.5 + quat1;
    check_iterable_approx<Quaternion<double>>::apply(
        quat3, Quaternion<double>(-8.7, -8.5, 10.9, -8.1));
    check_iterable_approx<Quaternion<double>>::apply(
        quat4, Quaternion<double>(-8.7, -8.5, 10.9, -8.1));
    // checks unary operator +
    check_iterable_approx<Quaternion<double>>::apply(quat3, +quat4);
    check_iterable_approx<Quaternion<double>>::apply(
        quat5, Quaternion<double>(3.6, 1.2, 1.3, 1.4));
    check_iterable_approx<Quaternion<double>>::apply(
        quat6, Quaternion<double>(3.6, 1.2, 1.3, 1.4));
    check_iterable_approx<Quaternion<double>>::apply(quat5, +quat6);

    quat3 = quat1 - quat2;
    quat4 = quat2 - quat1;
    quat5 = quat1 - 2.5;
    quat6 = 2.5 - quat1;
    check_iterable_approx<Quaternion<double>>::apply(
        quat3, Quaternion<double>(10.9, 10.9, -8.3, 10.9));
    check_iterable_approx<Quaternion<double>>::apply(
        quat4, Quaternion<double>(-10.9, -10.9, 8.3, -10.9));
    // checks unary operator -
    check_iterable_approx<Quaternion<double>>::apply(quat3, -quat4);
    check_iterable_approx<Quaternion<double>>::apply(
        quat5, Quaternion<double>(-1.4, 1.2, 1.3, 1.4));
    check_iterable_approx<Quaternion<double>>::apply(
        quat6, Quaternion<double>(1.4, -1.2, -1.3, -1.4));
    check_iterable_approx<Quaternion<double>>::apply(quat5, -quat6);

    quat3 = quat1 * quat2;
    quat4 = quat2 * quat1;
    quat5 = quat1 * 2.5;
    quat6 = 2.5 * quat1;
    check_iterable_approx<Quaternion<double>>::apply(
        quat3, Quaternion<double>(1.68, -48.22, -4.36, -0.04));
    check_iterable_approx<Quaternion<double>>::apply(
        quat4, Quaternion<double>(1.68, 3.36, 0.0, -48.3));
    CHECK(quat3 != quat4);
    check_iterable_approx<Quaternion<double>>::apply(
        quat5, Quaternion<double>(2.75, 3.0, 3.25, 3.5));
    check_iterable_approx<Quaternion<double>>::apply(
        quat6, Quaternion<double>(2.75, 3.0, 3.25, 3.5));
    CHECK(quat5 == quat6);

    quat3 = quat1 / quat2;
    quat4 = quat2 / quat1;
    quat5 = quat1 / 2.5;
    check_iterable_approx<Quaternion<double>>::apply(quat3,
                                                     quat1 * quat2.inv());
    check_iterable_approx<Quaternion<double>>::apply(quat4,
                                                     quat2 * quat1.inv());
    CHECK(quat3 != quat4);
    check_iterable_approx<Quaternion<double>>::apply(
        quat5, Quaternion<double>(0.44, 0.48, 0.52, 0.56));
  }

  {
    INFO("Check += -= *= /= operators")
    Quaternion<double> quat1(1.1, 1.2, 1.3, 1.4);
    Quaternion<double> quat1reset(1.1, 1.2, 1.3, 1.4);
    Quaternion<double> quat2(-9.8, -9.7, 9.6, -9.5);

    quat1 += quat2;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, Quaternion<double>(-8.7, -8.5, 10.9, -8.1));
    quat1 += 2.5;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, Quaternion<double>(-6.2, -8.5, 10.9, -8.1));
    quat1 = quat1reset;

    quat1 -= quat2;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, Quaternion<double>(10.9, 10.9, -8.3, 10.9));
    quat1 -= 2.5;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, Quaternion<double>(8.4, 10.9, -8.3, 10.9));
    quat1 = quat1reset;

    quat1 *= quat2;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, Quaternion<double>(1.68, -48.22, -4.36, -0.04));
    quat1 *= 2.5;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, Quaternion<double>(4.2, -120.55, -10.9, -0.1));
    quat1 = quat1reset;

    quat1 /= quat2;
    check_iterable_approx<Quaternion<double>>::apply(quat1,
                                                     quat1reset * quat2.inv());
    quat1 /= 2.5;
    check_iterable_approx<Quaternion<double>>::apply(
        quat1, quat1reset * quat2.inv() / 2.5);
  }
}
