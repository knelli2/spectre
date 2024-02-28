// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <variant>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/IdPair.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/CallWithDynamicType.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Frame {
struct Grid;
struct Distorted;
struct Inertial;
struct NoFrame;
struct BlockLogical;
}  // namespace Frame
/// \endcond

namespace intrp2::points {
/*!
 * \brief Base class for computing the target points we would like to
 * interpolate tensors to.
 */
template <size_t Dim>
class Target : public PUP::able {
 public:
  using BlockCoords = std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>;

  /// \cond
  Target() = default;
  Target(const Target&) = default;
  Target(Target&&) = default;
  Target& operator=(const Target&) = default;
  Target& operator=(Target&&) = default;
  /// \endcond

  ~Target() override = default;
  explicit Target(CkMigrateMessage* msg) : PUP::able(msg) {}
  WRAPPED_PUPable_abstract(Target);  // NOLINT

  /*!
   * \brief Name of this target as a string.
   */
  virtual const std::string& name() const = 0;

  using argument_tags = tmpl::list<>;

  template <typename Metavariables>
  std::optional<BlockCoords> operator()(
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& frame) const {
    return {block_logical_coordinates_in_frame(
        cache, time, this->target_points_no_frame(), frame)};
  }

  /*!
   * \brief Returns the block logical coordinates of the target points.
   *
   * In order to allow for a runtime choice of frame using this single base
   * class, a derived class must privately override the following virtual
   * function:
   *
   * \code {.cpp}
   * const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
   *     override;
   * \endcode
   *
   * This means that inside the derived class, the target points should not be
   * stored in a tensor of the proper frame (because the derived class won't
   * know the proper frame at compile time). But rather in a tensor of the
   * proper dimension and `Frame::NoFrame`. This then allows us to choose the
   * frame at runtime.
   */
  template <typename DbTags, typename Metavariables>
  std::optional<BlockCoords> block_logical_coordinates(
      const db::DataBox<DbTags>& box,
      const Parallel::GlobalCache<Metavariables>& cache, const double time,
      const std::string& frame) {
    using factory_classes =
        typename std::decay_t<Metavariables>::factory_creation::factory_classes;
    return call_with_dynamic_type<BlockCoords,
                                  tmpl::at<factory_classes, Target>>(
        this, [&](auto* const target) {
          return db::apply(*target, box, cache, time, frame);
        });
  }

  virtual const tnsr::I<DataVector, Dim, Frame::NoFrame>&
  target_points_no_frame() const = 0;
};
}  // namespace intrp2::points
