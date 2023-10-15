// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines MathFunctions::Sinusoid.

#pragma once

#include <memory>
#include <pup.h>

#include "Options/String.hpp"
#include "PointwiseFunctions/MathFunctions/MathFunction.hpp"  // IWYU pragma: keep
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

namespace MathFunctions {
template <size_t VolumeDim, typename Fr>
class Sinusoid;

/*!
 *  \ingroup MathFunctionsGroup
 *  \brief Sinusoid \f$f = A \sin\left(k x + \delta \right)\f$.
 *
 *  \details Input file options are: Amplitude, Phase, and Wavenumber
 */
template <typename Fr>
class Sinusoid<1, Fr> : public MathFunction<1, Fr> {
 public:
  struct Amplitude {
    using type = double;
    inline const static std::string help {"The amplitude."};
  };

  struct Wavenumber {
    using type = double;
    inline const static std::string help {"The wavenumber."};
  };

  struct Phase {
    using type = double;
    inline const static std::string help {"The phase shift."};
  };
  using options = tmpl::list<Amplitude, Wavenumber, Phase>;

  inline const static std::string help {
      "Applies a Sinusoid function to the input value"};

  Sinusoid(double amplitude, double wavenumber, double phase);
  Sinusoid() = default;
  std::unique_ptr<MathFunction<1, Fr>> get_clone() const override;

  WRAPPED_PUPable_decl_base_template(SINGLE_ARG(MathFunction<1, Fr>),
                                     Sinusoid);  // NOLINT

  explicit Sinusoid(CkMigrateMessage* /*unused*/) {}

  double operator()(const double& x) const override;
  DataVector operator()(const DataVector& x) const override;

  double first_deriv(const double& x) const override;
  DataVector first_deriv(const DataVector& x) const override;

  double second_deriv(const double& x) const override;
  DataVector second_deriv(const DataVector& x) const override;

  double third_deriv(const double& x) const override;
  DataVector third_deriv(const DataVector& x) const override;

  bool operator==(const MathFunction<1, Fr>& other) const override;
  bool operator!=(const MathFunction<1, Fr>& other) const override;
  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  double amplitude_{};
  double wavenumber_{};
  double phase_{};

  template <typename T>
  T apply_call_operator(const T& x) const;
  template <typename T>
  T apply_first_deriv(const T& x) const;
  template <typename T>
  T apply_second_deriv(const T& x) const;
  template <typename T>
  T apply_third_deriv(const T& x) const;
};

}  // namespace MathFunctions

/// \cond
template <typename Fr>
PUP::able::PUP_ID MathFunctions::Sinusoid<1, Fr>::my_PUP_ID = 0;  // NOLINT
/// \endcond
