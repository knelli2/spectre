// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

namespace intrp2::protocols {
/*!
 * \brief A protocol for the type alias `points` in a `intrp2::Target`.
 *
 * \details A struct conforming to the `Points` protocol must
 *
 * - have a type alias `tags_on_target` which is any simple or compute tags that
 *   these points offer in the target DataBox.
 *
 * - have a type alias `points_volume_compute_tags` that is a list of compute
 * tags in the element DataBox that will be used to compute the points.
 *
 * - be option creatable.
 *
 * - be pup-able.
 *
 * - have a function that returns the points with signature
 *
 * \code{cpp}
 * const tnsr::I<DataVector, Dim, Frame::NoFrame>&
 *     target_points_no_frame() const;
 * \endcode
 *
 * - have a function `size_t number_of_sets_of_points() const` that returns how
 *   many sets of points are stored in the `DataVector`s of the tensor returned
 *   from the `target_points_no_frame()`. This is useful, for example, if you
 *   want multiple concentric spheres in the same target.
 *
 * A struct that conforms to this protocol can also optionally conform to the
 * `db:protocols::Apply` protocol:
 *
 * TODO: Add example
 */
struct Points {
  template <typename ConformingType>
  struct test {
    using tags_on_target = typename ConformingType::tags_on_target;

    using points_volume_compute_tags =
        typename ConformingType::points_volume_compute_tags;

    using options = typename ConformingType::options;
  };
};
}  // namespace intrp2::protocols
