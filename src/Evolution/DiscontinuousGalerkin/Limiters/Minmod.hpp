// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "DataStructures/Tags.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/DiscontinuousGalerkin/Limiters/MinmodImpl.hpp"
#include "Evolution/DiscontinuousGalerkin/Limiters/MinmodType.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
class DataVector;
template <size_t VolumeDim>
class Direction;
template <size_t Dim, typename T>
class DirectionMap;
template <size_t VolumeDim>
class Element;
template <size_t VolumeDim>
class ElementId;
template <size_t VolumeDim>
class Mesh;
template <size_t VolumeDim>
class OrientationMap;

namespace boost {
template <class T>
struct hash;
}  // namespace boost

namespace PUP {
class er;
}  // namespace PUP

namespace Limiters::Minmod_detail {
template <size_t VolumeDim>
class BufferWrapper;
}  // namespace Limiters::Minmod_detail

namespace domain::Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
template <size_t VolumeDim>
struct Element;
template <size_t VolumeDim>
struct Mesh;
template <size_t VolumeDim>
struct SizeOfElement;
}  // namespace domain::Tags
/// \endcond

namespace Limiters {
/// \ingroup LimitersGroup
/// \brief A minmod-based generalized slope limiter
///
/// Implements the three minmod-based generalized slope limiters from
/// \cite Cockburn1999 Sec. 2.4: \f$\Lambda\Pi^1\f$, \f$\Lambda\Pi^N\f$, and
/// MUSCL. The implementation is system-agnostic and can act on an arbitrary
/// set of tensors.
///
/// #### Summary of the generalized slope limiter algorithms:
///
/// The MUSCL and \f$\Lambda\Pi^1\f$ limiters are both intended for use on
/// piecewise-linear solutions, i.e., on linear-order elements with two points
/// per dimension. These limiters operate by reducing the spatial slope of the
/// tensor components if the data look like they may contain oscillations.
/// Between these two, MUSCL is more dissipative --- it more aggressively
/// reduces the slopes of the data, so it may better handle strong shocks, but
/// it correspondingly produces the most broadening of features in the solution.
///
/// Note that we do not _require_ the MUSCL and \f$\Lambda\Pi^1\f$ limiters to
/// be used on linear-order elements. However, when they are used on a
/// higher-resolution grid, the limiters act to linearize the solution (by
/// discarding higher-order mode content) whether or not the slopes must be
/// reduced.
///
/// The \f$\Lambda\Pi^N\f$ limiter is intended for use with higher-order
/// elements (with more than two points per dimension), where the solution is a
/// piecewise polynomial of higher-than-linear order. This limiter generalizes
/// \f$\Lambda\Pi^1\f$: the post-limiter solution is the linearized solution of
/// \f$\Lambda\Pi^1\f$ in the case that the slopes must be reduced, but is the
/// original (higher-order) data in the case that the slopes are acceptable.
///
/// For all three types of minmod limiter, the algorithm can be relaxed from
/// TVD (total variation diminishing) in the means to TVB (total variation
/// bound) in the means. This may avoid limiting away smooth extrema in the
/// solution that would otherwise look like spurious oscillations. When this
/// correction is enabled, the limiter will not reduce the slope (but may still
/// linearize) on elements where the slope is less than \f$m h^2\f$, where
/// \f$m\f$ is the TVB constant and \f$h\f$ is the size of the DG element.
/// Note the "in the means" qualifier: the limiter controls the oscillation
/// between the mean solution values across neighboring cells, but may not
/// control oscillations within the cells.
///
/// The choice of the TVB constant \f$m\f$ is difficult. Larger values result in
/// fewer limiter activations, especially near smooth extrema in the solution
/// --- this can help to avoid incorrectly limiting away these smooth extrema,
/// but can also result in insufficient limiting of truly spurious oscillations.
/// The reference uses a value of 50 when presenting the limiter with simple
/// shock tests, but in general the value \f$m\f$ that optimizes between
/// robustness and minimal loss of accuracy is problem dependent.
///
/// #### Notes on the SpECTRE implementation of the generalized slope limiters:
///
/// This implementation can act on an arbitrary set of tensors; the limiting
/// algorithm is applied to each component of each tensor independently. This
/// is a convenient and general interface. However, when the evolution system
/// has multiple evolved variables, the recommendation of the reference is to
/// apply the limiter to the system's characteristic variables to reduce
/// spurious post-limiting oscillations. In SpECTRE, applying the limiter to
/// the characteristic variables requires specializing the limiter to each
/// evolution system.
///
/// The limiter acts in the `Frame::ElementLogical` coordinates, because in
/// these coordinates it is straightforward to formulate the algorithm. This
/// means the limiter can operate on generic deformed grids --- however, some
/// things can start to break down, especially on strongly deformed grids:
/// 1. When the Jacobian (from `Frame::ElementLogical` to `Frame::Inertial`)
///    varies across the element, then the limiter fails to be conservative.
///    This is because the integral of a tensor `u` over the element will change
///    after the limiter activates on `u`.
/// 2. When there is a sudden change in the size of the elements (perhaps at an
///    h-refinement boundary, or at the boundary between two blocks with very
///    different mappings), a smooth solution in `Frame::Inertial` can appear
///    to have a kink in `Frame::ElementLogical`. The Minmod implementation
///    includes some (tested but unproven) corrections based on the size of the
///    elements that try to reduce spurious limiter activations near these fake
///    kinks.
///
/// When an element has multiple neighbors in any direction, an effective mean
/// and neighbor size in this direction are computed by averaging over the
/// multiple neighbors. This simple generalization of the minmod limiter enables
/// it to operate on h-refined grids.
///
/// \tparam VolumeDim The number of spatial dimensions.
/// \tparam TagsToLimit A typelist of tags specifying the tensors to limit.
template <size_t VolumeDim, typename TagsToLimit>
class Minmod;

template <size_t VolumeDim, typename... Tags>
class Minmod<VolumeDim, tmpl::list<Tags...>> {
 public:
  /// \brief The MinmodType
  ///
  /// One of `Limiters::MinmodType`. See `Limiters::Minmod`
  /// documentation for details. Note in particular that on grids with more than
  /// two points per dimension, the recommended type is `LambdaPiN`.
  struct Type {
    using type = MinmodType;
    inline const static std::string help {"Type of minmod"};
  };
  /// \brief The TVB constant
  ///
  /// See `Limiters::Minmod` documentation for details. The optimal value of
  /// this parameter is unfortunately problem-dependent.
  struct TvbConstant {
    using type = double;
    static type lower_bound() { return 0.0; }
    inline const static std::string help {"TVB constant 'm'"};
  };
  /// \brief Turn the limiter off
  ///
  /// This option exists to temporarily disable the limiter for debugging
  /// purposes. For problems where limiting is not needed, the preferred
  /// approach is to not compile the limiter into the executable.
  struct DisableForDebugging {
    using type = bool;
    static type suggested_value() { return false; }
    inline const static std::string help {"Disable the limiter"};
  };
  using options = tmpl::list<Type, TvbConstant, DisableForDebugging>;
  inline const static std::string help {
      "A minmod-based generalized slope limiter.\n"
      "The different types of minmod are more or less aggressive in trying\n"
      "to reduce slopes. The TVB correction allows the limiter to ignore\n"
      "'small' slopes, and helps to avoid limiting of smooth extrema in the\n"
      "solution.\n"};

  /// \brief Constuct a Minmod slope limiter
  ///
  /// \param minmod_type The type of Minmod slope limiter.
  /// \param tvb_constant The value of the TVB constant.
  /// \param disable_for_debugging Switch to turn the limiter off.
  explicit Minmod(MinmodType minmod_type, double tvb_constant,
                  bool disable_for_debugging = false);

  Minmod() = default;
  Minmod(const Minmod& /*rhs*/) = default;
  Minmod& operator=(const Minmod& /*rhs*/) = default;
  Minmod(Minmod&& /*rhs*/) = default;
  Minmod& operator=(Minmod&& /*rhs*/) = default;
  ~Minmod() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p);

  // To facilitate testing
  /// \cond
  MinmodType minmod_type() const { return minmod_type_; }
  /// \endcond

  /// \brief Data to send to neighbor elements.
  struct PackagedData {
    tuples::TaggedTuple<::Tags::Mean<Tags>...> means;
    std::array<double, VolumeDim> element_size =
        make_array<VolumeDim>(std::numeric_limits<double>::signaling_NaN());

    // NOLINTNEXTLINE(google-runtime-references)
    void pup(PUP::er& p) {
      p | means;
      p | element_size;
    }
  };

  using package_argument_tags =
      tmpl::list<Tags..., domain::Tags::Mesh<VolumeDim>,
                 domain::Tags::SizeOfElement<VolumeDim>>;

  /// \brief Package data for sending to neighbor elements.
  ///
  /// The following quantities are stored in `PackagedData` and communicated
  /// between neighboring elements:
  /// - the cell-averaged mean of each tensor component, and
  /// - the size of the cell along each logical coordinate direction.
  ///
  /// \param packaged_data The data package to fill with this element's values.
  /// \param tensors The tensors to be averaged and packaged.
  /// \param mesh The mesh on which the tensor values are measured.
  /// \param element_size The size of the element in inertial coordinates, along
  ///        each dimension of logical coordinates.
  /// \param orientation_map The orientation of the neighbor
  void package_data(gsl::not_null<PackagedData*> packaged_data,
                    const typename Tags::type&... tensors,
                    const Mesh<VolumeDim>& mesh,
                    const std::array<double, VolumeDim>& element_size,
                    const OrientationMap<VolumeDim>& orientation_map) const;

  using limit_tags = tmpl::list<Tags...>;
  using limit_argument_tags =
      tmpl::list<domain::Tags::Mesh<VolumeDim>,
                 domain::Tags::Element<VolumeDim>,
                 domain::Tags::Coordinates<VolumeDim, Frame::ElementLogical>,
                 domain::Tags::SizeOfElement<VolumeDim>>;

  /// \brief Limits the solution on the element.
  ///
  /// For each component of each tensor, the limiter will (in general) linearize
  /// the data, then possibly reduce its slope, dimension-by-dimension, until it
  /// no longer looks oscillatory.
  ///
  /// \param tensors The tensors to be limited.
  /// \param mesh The mesh on which the tensor values are measured.
  /// \param element The element on which the tensors to limit live.
  /// \param logical_coords The element logical coordinates of the mesh
  ///        gridpoints.
  /// \param element_size The size of the element, in the inertial coordinates.
  /// \param neighbor_data The data from each neighbor.
  ///
  /// \return whether the limiter modified the solution or not.
  ///
  /// \note The return value is false if the limiter knows it has not modified
  /// the solution. True return values can indicate:
  /// - The solution was limited to reduce the slope, whether by a large factor
  ///   or by a factor only roundoff away from unity.
  /// - The solution was linearized but not limited.
  /// - The solution is identical to the input, if the input was a linear
  ///   function on a higher-order mesh, so that the limiter cannot know that
  ///   the linearization step did not actually modify the data. This is
  ///   somewhat contrived and is unlikely to occur outside of code tests or
  ///   test cases with very clean initial data.
  bool operator()(
      const gsl::not_null<std::add_pointer_t<typename Tags::type>>... tensors,
      const Mesh<VolumeDim>& mesh, const Element<VolumeDim>& element,
      const tnsr::I<DataVector, VolumeDim, Frame::ElementLogical>&
          logical_coords,
      const std::array<double, VolumeDim>& element_size,
      const std::unordered_map<
          std::pair<Direction<VolumeDim>, ElementId<VolumeDim>>, PackagedData,
          boost::hash<std::pair<Direction<VolumeDim>, ElementId<VolumeDim>>>>&
          neighbor_data) const;

 private:
  template <size_t LocalDim, typename LocalTagList>
  // NOLINTNEXTLINE(readability-redundant-declaration) false positive
  friend bool operator==(const Minmod<LocalDim, LocalTagList>& lhs,
                         const Minmod<LocalDim, LocalTagList>& rhs);

  MinmodType minmod_type_;
  double tvb_constant_;
  bool disable_for_debugging_;
};

template <size_t VolumeDim, typename... Tags>
Minmod<VolumeDim, tmpl::list<Tags...>>::Minmod(const MinmodType minmod_type,
                                               const double tvb_constant,
                                               const bool disable_for_debugging)
    : minmod_type_(minmod_type),
      tvb_constant_(tvb_constant),
      disable_for_debugging_(disable_for_debugging) {
  ASSERT(tvb_constant >= 0.0, "The TVB constant must be non-negative.");
}

template <size_t VolumeDim, typename... Tags>
void Minmod<VolumeDim, tmpl::list<Tags...>>::pup(PUP::er& p) {
  p | minmod_type_;
  p | tvb_constant_;
  p | disable_for_debugging_;
}

template <size_t VolumeDim, typename... Tags>
void Minmod<VolumeDim, tmpl::list<Tags...>>::package_data(
    const gsl::not_null<PackagedData*> packaged_data,
    const typename Tags::type&... tensors, const Mesh<VolumeDim>& mesh,
    const std::array<double, VolumeDim>& element_size,
    const OrientationMap<VolumeDim>& orientation_map) const {
  if (UNLIKELY(disable_for_debugging_)) {
    // Do not initialize packaged_data
    return;
  }

  const auto wrap_compute_means = [&mesh, &packaged_data](auto tag,
                                                          const auto tensor) {
    for (size_t i = 0; i < tensor.size(); ++i) {
      // Compute the mean using the local orientation of the tensor and mesh:
      // this avoids the work of reorienting the tensor while giving the same
      // result.
      get<::Tags::Mean<decltype(tag)>>(packaged_data->means)[i] =
          mean_value(tensor[i], mesh);
    }
    return '0';
  };
  expand_pack(wrap_compute_means(Tags{}, tensors)...);
  packaged_data->element_size =
      orientation_map.permute_from_neighbor(element_size);
}

template <size_t VolumeDim, typename... Tags>
bool Minmod<VolumeDim, tmpl::list<Tags...>>::operator()(
    const gsl::not_null<std::add_pointer_t<typename Tags::type>>... tensors,
    const Mesh<VolumeDim>& mesh, const Element<VolumeDim>& element,
    const tnsr::I<DataVector, VolumeDim, Frame::ElementLogical>& logical_coords,
    const std::array<double, VolumeDim>& element_size,
    const std::unordered_map<
        std::pair<Direction<VolumeDim>, ElementId<VolumeDim>>, PackagedData,
        boost::hash<std::pair<Direction<VolumeDim>, ElementId<VolumeDim>>>>&
        neighbor_data) const {
  if (UNLIKELY(disable_for_debugging_)) {
    // Do not modify input tensors
    return false;
  }

  DataVector u_lin_buffer(mesh.number_of_grid_points());
  Minmod_detail::BufferWrapper<VolumeDim> buffer(mesh);

  bool limiter_activated = false;
  const auto wrap_minmod_impl = [this, &limiter_activated, &element, &mesh,
                                 &logical_coords, &element_size, &neighbor_data,
                                 &u_lin_buffer,
                                 &buffer](auto tag, const auto tensor) {
    limiter_activated =
        Minmod_detail::minmod_impl<VolumeDim, decltype(tag)>(
            &u_lin_buffer, &buffer, tensor, minmod_type_, tvb_constant_, mesh,
            element, logical_coords, element_size, neighbor_data) or
        limiter_activated;
    return '0';
  };
  expand_pack(wrap_minmod_impl(Tags{}, tensors)...);
  return limiter_activated;
}

template <size_t LocalDim, typename LocalTagList>
bool operator==(const Minmod<LocalDim, LocalTagList>& lhs,
                const Minmod<LocalDim, LocalTagList>& rhs) {
  return lhs.minmod_type_ == rhs.minmod_type_ and
         lhs.tvb_constant_ == rhs.tvb_constant_ and
         lhs.disable_for_debugging_ == rhs.disable_for_debugging_;
}

template <size_t VolumeDim, typename TagList>
bool operator!=(const Minmod<VolumeDim, TagList>& lhs,
                const Minmod<VolumeDim, TagList>& rhs) {
  return not(lhs == rhs);
}

}  // namespace Limiters
