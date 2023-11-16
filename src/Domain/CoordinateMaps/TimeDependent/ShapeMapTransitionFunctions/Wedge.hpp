// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/Falloff.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/Structure/Direction.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {

/*!
 * \brief A transition function that falls off from 1 at an inner surface of
 * a wedge to 0 at an outer surface of a wedge. Meant to be used in
 * `domain::CoordinateMaps::Wedge` blocks.
 *
 * \details There are two functional forms that this transition can take
 * depending on the \p falloff parameter passed to the constructor. If the
 * falloff is
 * `domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff::Linear` then
 * the functional form of this transition is
 *
 * \begin{equation}
 * f(r, \theta, \phi) = \frac{D_{\text{out}}(r, \theta, \phi) -
 * r}{D_{\text{out}}(r, \theta, \phi) - D_{\text{in}}(r, \theta, \phi)}.
 * \label{eq:linear_transition}
 * \end{equation}
 *
 * If the falloff is
 * `domain::CoordinateMaps::ShapeMapTransitionFunctions::Falloff::Inverse`, then
 * the functional form is
 *
 * \begin{equation}
 * f(r, \theta, \phi) = a + \frac{b}{r},
 * \label{eq:inverse_transition}
 * \end{equation}
 *
 * where
 *
 * \f{align}{
 * a &= \frac{-D_{\text{in}}(r, \theta, \phi)}{D_{\text{out}}(r, \theta, \phi) -
 * D_{\text{in}}(r, \theta, \phi)} \\
 * b &= \frac{D_{\text{out}}(r, \theta, \phi)D_{\text{in}}(r, \theta,
 * \phi)}{D_{\text{out}}(r, \theta, \phi) -
 * D_{\text{in}}(r, \theta, \phi)} = -aD_{\text{out}}(r, \theta, \phi).
 * \label{eq:inverse_a_and_b}
 * \f}
 *
 * For both of these falloffs, we have that
 *
 * \begin{equation}
 * D(r, \theta, \phi) = R\left(\frac{1 -
 * s}{\sqrt{3}\max(|\sin{\theta}\cos{\phi}|,|\sin{\theta}\sin{\phi}|,
 * |\cos{\theta}|)} + s\right),
 * \label{eq:distance}
 * \end{equation}
 *
 * where $s$ is the sphericity of the surface which goes from 0 (flat) to 1
 * (spherical), $R$ is the radius of the spherical surface, $\text{out}$ is the
 * outer surface, and $\text{in}$ is the inner surface. If the sphericity is 1,
 * then the surface is a sphere so $D = R$. If the sphericity is 0, then the
 * surface is a cube. This cube is circumscribed by a sphere of radius $R$. See
 * `domain::CoordinateMaps::Wedge` for more of an explanation of these boundary
 * surfaces and their sphericities.
 *
 * \note Because the shape map distorts only radii and does not affect angles,
 * $D$ is only a function of $\theta$ so we have that $D(r,\theta,\phi) =
 * D(\theta)$.
 *
 * There are several assumptions made for this mapping:
 *
 * - The `domain::CoordinateMaps::Wedge`s have an opening angle of
 *   $\frac{\pi}{2}$ in each direction.
 * - The coordinates $r, \theta, \phi$ are assumed to be from the center of the
 *   wedge, not the center of the computational domain.
 * - The wedges are concentric. (see the constructor)
 * - The $\max$ in the denominator of $\ref{eq:distance}$ can be simplified a
 *   bit to $\max(|x|, |y|, |z|)/r$. It was written the other way in
 *   $\ref{eq:distance}$ to emphasize that $D$ has no radial dependence.
 *
 * ## Gradient
 *
 * The cartesian gradient of the linear transition function (Eq.
 * $\ref{eq:linear_transition}$) is
 *
 * \begin{equation}
 * \frac{\partial f}{\partial x_i} = \frac{\frac{\partial
 * D_{\text{out}}}{\partial x_i} - \frac{x_i}{r}}{D_{\text{out}} -
 * D_{\text{in}}} - \frac{\left(D_{\text{out}} - r\right)\left(\frac{\partial
 * D_{\text{out}}}{\partial x_i} - \frac{\partial
 * D_{\text{in}}}{\partial x_i} \right)}{\left(D_{\text{out}} -
 * D_{\text{in}}\right)^2}.
 * \end{equation}
 *
 * The cartesian gradient of the inverse transition function (Eq.
 * $\ref{eq:inverse_transition}$) is
 *
 * \begin{equation}
 * \frac{\partial f}{\partial x_i} = \frac{\partial a}{\partial x_i} +
 * \frac{1}{r}\frac{\partial b}{\partial x_i} - \frac{bx_i}{r^3}
 * \end{equation}
 *
 * where
 *
 * \f{align}{
 * \frac{\partial a}{\partial x_i} &= \frac{D_{\text{in}}\frac{\partial
 * D_{\text{out}}}{\partial x_i} -D_{\text{out}}\frac{\partial
 * D_{\text{in}}}{\partial x_i} }{(D_{\text{out}} - D_{\text{in}})^2} \\
 * \frac{\partial b}{\partial x_i} &= -a\frac{\partial D_{\text{out}}}{\partial
 * x_i} - D_{\text{out}}\frac{\partial a}{\partial x_i}.
 * \f}
 *
 * The computation of the gradient of $D$ depends on the result of the $\max$ in
 * the denominator of $\ref{eq:distance}$. To simplify the expression, let $j
 * \in \{0,1,2\}$ correspond to the $\max(|x|, |y|, |z|)$ such that $j=0$ if
 * $|x|$ is the $\max$, and so on. Then we can write the gradient of $D$ as
 *
 * \f{align}{
 * \frac{\partial D}{\partial x_j} &= -\text{sgn}(x_j)\frac{R(1-s)\left(r^2 -
 * x_j^2\right)}{r x_j^2 \sqrt{3}} \\
 * \frac{\partial D}{\partial x_{j+1}} &= \frac{R(1-s)x_{j+1}}{r |x_j| \sqrt{3}}
 * \\
 * \frac{\partial D}{\partial x_{j+2}} &= \frac{R(1-s)x_{j+2}}{r |x_j| \sqrt{3}}
 * \f}
 *
 * where $j+1$ and $j+2$ are understood to be $\mod 3$. Also since we don't
 * allow points at the origin of these wedges, we can be assured that for
 * whichever $j$ is max, that $x_j$ won't be zero.
 *
 * ## Original radius divided by mapped radius
 *
 * Given an already mapped point and the distorted radius $\Sigma(\theta, \phi)
 * = \sum_{\ell,m}\lambda_{\ell,m}Y_{\ell,m}(\theta, \phi)$, we can figure out
 * the ratio of the original radius $r$ to the mapped $\tilde{r}$ by solving
 *
 * \begin{equation}
 * \frac{r}{\tilde{r}} = \frac{1}{1-f(r,\theta,\phi)\Sigma(\theta,\phi)}.
 * \end{equation}
 *
 * For the linear transition, this will result in us having to solve a quadratic
 * equation of the form
 *
 * \begin{equation}
 * \tilde{r} x^2 + \left(\frac{D_{\text{out}} -
 * D_{\text{in}}}{\Sigma(\theta,\phi)} - D_{\text{out}}\right) x -
 * \frac{D_{\text{out}} - D_{\text{in}}}{\Sigma(\theta,\phi)} = 0
 * \label{eq:linear_r_over_rtil}
 * \end{equation}
 *
 * where $x = \frac{r}{\tilde{r}}$.
 *
 * For the inverse transition, we get an analytic expression for the original
 * radius over the mapped radius of the form
 *
 * \begin{equation}
 * \frac{r}{\tilde{r}} = \frac{1 + \frac{b\Sigma(\theta,\phi)}{\tilde{r}}}{1 -
 * a\Sigma(\theta,\phi)}.
 * \label{eq:inverse_r_over_rtil}
 * \end{equation}
 *
 * \note Since $D$ is not a function of the radius, we can treat it as a
 * constant in Eqs. $\ref{eq:linear_r_over_rtil}$ and
 * $\ref{eq:inverse_r_over_rtil}$ and compute the angles using $\tilde{r}$.
 *
 * ## Special cases
 *
 * The equations above become simpler if the inner boundary, outer boundary, or
 * both are spherical ($s = 1$). If a boundary is spherical, then $D = R$ and
 * $\frac{\partial D}{\partial x_i} = 0$. This simplifies the expanded form of
 * the gradient of both falloffs significantly. Also, for computing
 * $\frac{r}{\tilde{r}$, since $D$ doesn't depend on angles we can use
 * $\tilde{r}$ to compute $D$ instead of $r$.
 */
class Wedge final : public ShapeMapTransitionFunction {
  struct Surface {
    double radius{};
    double sphericity{};

    // This is the distance from the center (assumed to be 0,0,0) to this
    // surface in the same direction as coords. The calculation is cheaper if
    // you know the axis ahead of time
    template <typename T>
    T distance(const std::array<T, 3>& coords,
               const std::optional<size_t>& axis = std::nullopt) const;

    void pup(PUP::er& p);

    bool operator==(const Surface& other) const;
    bool operator!=(const Surface& other) const;
  };

 public:
  explicit Wedge() = default;

  /*!
   * \brief Construct a Wedge transition for a wedge block in a given direction.
   *
   * \details Many concentric wedges can be part of the same falloff from 1 at
   * the inner boundary of the innermost wedge, to 0 at the outer boundary of
   * the outermost wedge.
   *
   * \param inner_radius Inner radius of innermost wedge
   * \param outer_radius Outermost radius of outermost wedge
   * \param inner_sphericity Sphericity of innermost surface of innermost wedge
   * \param outer_sphericity Sphericity of outermost surface of outermost wedge
   * \param falloff How the transition function falls off to zero
   * \param axis The direction that this wedge is in. Both the positive and
   * negative direction get the same axis.
   */
  Wedge(double inner_radius, double outer_radius, double inner_sphericity,
        double outer_sphericity, Falloff falloff, size_t axis);

  double operator()(const std::array<double, 3>& source_coords) const override;
  DataVector operator()(
      const std::array<DataVector, 3>& source_coords) const override;

  std::optional<double> original_radius_over_radius(
      const std::array<double, 3>& target_coords,
      double distorted_radius) const override;

  double map_over_radius(
      const std::array<double, 3>& source_coords) const override;
  DataVector map_over_radius(
      const std::array<DataVector, 3>& source_coords) const override;

  std::array<double, 3> gradient(
      const std::array<double, 3>& source_coords) const override;
  std::array<DataVector, 3> gradient(
      const std::array<DataVector, 3>& source_coords) const override;

  WRAPPED_PUPable_decl_template(Wedge);
  explicit Wedge(CkMigrateMessage* msg);
  void pup(PUP::er& p) override;

  std::unique_ptr<ShapeMapTransitionFunction> get_clone() const override {
    return std::make_unique<Wedge>(*this);
  }

  bool operator==(const ShapeMapTransitionFunction& other) const override;
  bool operator!=(const ShapeMapTransitionFunction& other) const override;

 private:
  template <typename T>
  T call_impl(const std::array<T, 3>& source_coords) const;

  template <typename T>
  T map_over_radius_impl(const std::array<T, 3>& source_coords) const;

  template <typename T>
  std::array<T, 3> gradient_impl(const std::array<T, 3>& source_coords) const;

  template <typename T>
  void check_distances(const std::array<T, 3>& coords) const;

  Surface inner_surface_{};
  Surface outer_surface_{};
  Falloff falloff_{};
  size_t axis_{};
  static constexpr double eps_ = std::numeric_limits<double>::epsilon() * 100;
};
}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
