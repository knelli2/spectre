// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/SpinWeighted.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/Cce/AnalyticSolutions/WorldtubeData.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

namespace Cce::Solutions {

/*!
 * \brief Analytic solution representing a coordinate oscillation about a
 * stationary Schwarzschild black hole.
 *
 * \details As the oscillation in the metric data at the worldtube is a pure
 * coordinate effect, the system evolved using this worldtube data should
 * produce zero news. The solution is a coordinate transform applied to the
 * Schwarzschild solution in Kerr-Schild coordinates.
 */
struct BouncingBlackHole : public WorldtubeData {
  struct Amplitude {
    using type = double;
    inline const static std::string help{
        "The coordinate distance of the gauge oscillation"};
    static type lower_bound() { return 0.0; }
    static type suggested_value() { return 2.0; }
  };
  struct ExtractionRadius {
    using type = double;
    inline const static std::string help{
        "The extraction radius of the spherical solution"};
    static type lower_bound() { return 0.0; }
    static type suggested_value() { return 20.0; }
  };
  struct Mass {
    using type = double;
    inline const static std::string help{
        "The mass of the Schwarzschild black hole"};
    static type lower_bound() { return 0.0; }
    static type suggested_value() { return 1.0; }
  };
  struct Period {
    using type = double;
    inline const static std::string help{
        "The period of the coordinate oscillation"};
    static type lower_bound() { return 0.0; }
    static type suggested_value() { return 40.0; }
  };

  inline const static std::string help{
      "Analytic solution in which a static black hole is placed in an "
      "oscillating coordinate system"};

  using options = tmpl::list<Amplitude, ExtractionRadius, Mass, Period>;

  WRAPPED_PUPable_decl_template(BouncingBlackHole);  // NOLINT

  explicit BouncingBlackHole(CkMigrateMessage* msg) : WorldtubeData(msg) {}

  // clang doesn't manage to use = default correctly in this case
  // NOLINTNEXTLINE(modernize-use-equals-default)
  BouncingBlackHole() {}

  BouncingBlackHole(double amplitude, double extraction_radius, double mass,
                    double period);

  std::unique_ptr<WorldtubeData> get_clone() const override;

  void pup(PUP::er& p) override;

 protected:
  // The bouncing black hole solution is easily computed directly, so requires
  // no additional preparation.
  void prepare_solution(const size_t /*l_max*/,
                        const double /*time*/) const override{};

  using WorldtubeData::variables_impl;

  /*!
   * \brief The implementation function that computes the spacetime metric on
   * the extraction sphere at collocation points associated with angular
   * resolution `l_max`.
   *
   * \details The spacetime metric \f$g_{a b}\f$ is determined by evaluating the
   * Kerr-Schild metric at a set of transformed coordinates \f$t^\prime = t,
   * y^\prime = y, z^\prime = z\f$, and
   *
   * \f{align*}{
   * x = x^\prime + A \left(\sin\left(\frac{2 \pi t}{T}\right)\right)^4,
   * \f}
   *
   * where the amplitude \f$A\f$ is set by the option `Amplitude` and the period
   * \f$T\f$ is set by the option `Period`. In this notation we take
   * the primed coordinates to be the coordinates for which the black hole has
   * time-dependent coordinate position.
   */
  void variables_impl(
      gsl::not_null<tnsr::aa<DataVector, 3>*> spacetime_metric, size_t l_max,
      double time,
      tmpl::type_<gr::Tags::SpacetimeMetric<DataVector, 3>> /*meta*/)
      const override;

  /*!
   * \brief The implementation function that computes the first time derivative
   * of the spacetime metric on the extraction sphere.
   *
   * \details The time derivative of the spacetime metric
   * \f$\partial_t g_{a b}\f$ comes entirely from the Jacobian factor:
   *
   * \f{align*}{
   * \partial_t x = \frac{8 \pi A}{T} \cos\left(\frac{2 \pi t}{T}\right)
   * \left(\sin\left(\frac{2 \pi t}{T}\right)\right)^3,
   * \f}
   *
   * so the transformed metric derivative is,
   *
   * \f{align*}{
   * \partial_t g_{a^\prime b^\prime} = 2 \partial_{(a^\prime} \partial_t x
   * \partial_{b^\prime)} x^a g_{x a}.
   * \f}
   *
   * In this notation we take the primed coordinates to be the coordinates for
   * which the black hole has time-dependent coordinate position.
   */
  void variables_impl(
      gsl::not_null<tnsr::aa<DataVector, 3>*> dt_spacetime_metric, size_t l_max,
      double time,
      tmpl::type_<
          ::Tags::dt<gr::Tags::SpacetimeMetric<DataVector, 3>>> /*meta*/)
      const override;

  /*!
   * \brief The implementation function that computes the first spatial
   * derivative of the spacetime metric on the extraction sphere.
   *
   * \details The calculation proceeds by standard coordinate transform
   * techniques for the transformation given by \f$t^\prime = t,
   * y^\prime = y, z^\prime = z\f$, and
   *
   * \f{align*}{
   * x = x^\prime + A \left(\sin\left(\frac{2 \pi t}{T}\right)\right)^4,
   * \f}
   *
   * The general coordinate transformation formula that gives the metric
   * is then
   * \f{align*}{
   * \partial_a g_{b c} =
   * \partial_a \partial_b x^{\prime a^\prime} \partial_c x^{\prime b^\prime}
   * g_{a^\prime b^\prime}
   * + \partial_b x^{\prime a^\prime} \partial_a \partial_c x^{\prime b^\prime}
   * g_{a^\prime b^\prime}
   * + \partial_a x^{\prime a^\prime} \partial_b x^{\prime b^\prime}
   * \partial_c x^{\prime c^\prime} \partial_a g_{b c}
   * \f}
   */
  void variables_impl(
      gsl::not_null<tnsr::iaa<DataVector, 3>*> d_spacetime_metric, size_t l_max,
      double time,
      tmpl::type_<gh::Tags::Phi<DataVector, 3>> /*meta*/) const override;

  /// The News in the bouncing black hole solution vanishes, as the oscillation
  /// comes entirely from a coordinate transform.
  void variables_impl(
      gsl::not_null<Scalar<SpinWeighted<ComplexDataVector, -2>>*> news,
      size_t output_l_max, double time,
      tmpl::type_<Tags::News> /*meta*/) const override;

  double amplitude_ = std::numeric_limits<double>::signaling_NaN();
  double mass_ = std::numeric_limits<double>::signaling_NaN();
  double frequency_ = std::numeric_limits<double>::signaling_NaN();
};
}  // namespace Cce::Solutions
