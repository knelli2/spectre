// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

namespace intrp2::protocols {
/*!
 * \brief A protocol for the callback done after interpolation has finished to a
 * target.
 *
 * \details A struct conforming to the `TargetPoints` protocol must have
 *
 * - a type alias `tags_to_observe_on_target` which is a list of simple or
 *   compute tags on the target corresponding to quantities that a user can
 *   specify in an input file to observe.
 *
 * - a type alias `non_observation_tags_on_target` which is a list of simple or
 *   compute tags on the target that are necessary to compute the
 *   `tags_to_observe_on_target`, but are not themselves observable. An example
 *   would be if you wanted to  observe the integral of a scalar over the
 *   points, then `non_observation_tags_on_target` would contain the compute tag
 *   for the scalar, but `tags_to_observe_on_target` would contain the compute
 *   tag for the integral of the scalar.
 *
 * - a type alias `volume_compute_tags` which is a list of compute tags in the
 *   element ObservationBox that will compute quantities before they are
 *   interpolated onto the target points.
 *
 * - inherit from the `intrp2::callbacks::Callback` abstract base class
 *
 * - be option creatable.
 *
 * - have a function `apply` with the following signature
 *
 * \code{cpp}
 * template <typename Metavariables, typename TemporalId>
 * static void apply(const db::Access& access,
 *                   Parallel::GlobalCache<Metavariables>& cache,
 *                   const TemporalId& time);
 * \endcode
 *
 * TODO: Add example
 */
struct Callback {
  template <typename ConformingType>
  struct test {
    using tags_to_observe_on_target =
        typename ConformingType::tags_to_observe_on_target;
    using non_observation_tags_on_target =
        typename ConformingType::non_observation_tags_on_target;

    using volume_compute_tags = typename ConformingType::volume_compute_tags;

    using options = typename ConformingType::options;
  };
};
}  // namespace intrp2::protocols
