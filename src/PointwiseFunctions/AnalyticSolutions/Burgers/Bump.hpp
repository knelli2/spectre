// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>

#include "DataStructures/DataBox/Prefixes.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/Burgers/Tags.hpp"  // IWYU pragma: keep
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
class DataVector;
// IWYU pragma: no_forward_declare Tensor
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace Burgers::Solutions {
/*!
 * \brief A solution resembling a bump.
 *
 * At \f$t=0\f$, the solution is a parabola:
 * \f{equation*}
 *  u(x, t) = h \left(1 - \left(\frac{x - c}{w}\right)^2\right),
 * \f}
 * where \f$h\f$ is the height, \f$c\f$ is the center, and \f$w\f$ is
 * the distance from the center to the zeros.  A shock propagates in
 * from infinity and reaches one of the zeros at \f$t = \frac{w}{2
 * h}\f$.
 */
class Bump : public evolution::initial_data::InitialData,
             public MarkAsAnalyticSolution {
 public:
  struct HalfWidth {
    using type = double;
    inline const static std::string help{
        "The distance from the center to the zero of the bump"};
    static type lower_bound() { return 0.; }
  };

  struct Height {
    using type = double;
    inline const static std::string help{"The height of the bump"};
  };

  struct Center {
    using type = double;
    inline const static std::string help{"The center of the bump"};
  };

  using options = tmpl::list<HalfWidth, Height, Center>;
  inline const static std::string help{"A bump solution"};

  Bump() = default;
  Bump(const Bump&) = default;
  Bump& operator=(const Bump&) = default;
  Bump(Bump&&) = default;
  Bump& operator=(Bump&&) = default;
  ~Bump() override = default;

  Bump(double half_width, double height, double center = 0.);

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit Bump(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Bump);
  /// \endcond

  template <typename T>
  Scalar<T> u(const tnsr::I<T, 1>& x, double t) const;

  template <typename T>
  Scalar<T> du_dt(const tnsr::I<T, 1>& x, double t) const;

  tuples::TaggedTuple<Tags::U> variables(const tnsr::I<DataVector, 1>& x,
                                         double t,
                                         tmpl::list<Tags::U> /*meta*/) const;

  tuples::TaggedTuple<::Tags::dt<Burgers::Tags::U>> variables(
      const tnsr::I<DataVector, 1>& x, double t,
      tmpl::list<::Tags::dt<Tags::U>> /*meta*/) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  double half_width_ = std::numeric_limits<double>::signaling_NaN();
  double height_ = std::numeric_limits<double>::signaling_NaN();
  double center_ = std::numeric_limits<double>::signaling_NaN();
};
}  // namespace Burgers::Solutions
