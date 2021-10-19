// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/OptionTags.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Options/Options.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace domain {
namespace creators {
namespace time_dependence {
/*!
 * \brief A TimeDependence that is a composition of two `UniformTranslation`s.
 *
 * To create this from options, use something like
 *
 * ```
 * CompositionUniformTranslation:
 *   UniformTranslation:
 *     UniformTranslation:
 *       InitialTime:
 *       Velocity:
 *   UniformTranslation1:
 *     UniformTranslation1:
 *       InitialTime:
 *       Velocity:
 * ```
 *
 * The reason for the `UniformTranslation1` is so that the factory can
 * distinguish the two and so that the names of the FunctionsOfTime are unique.
 */
template <size_t MeshDim>
class CompositionUniformTranslation final : public TimeDependence<MeshDim> {
 private:
  using TranslationMap =
      domain::CoordinateMaps::TimeDependent::Translation<MeshDim>;

 public:
  static constexpr size_t mesh_dim = MeshDim;

  using CoordMap = domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                         TranslationMap, TranslationMap>;

  using maps_list = tmpl::list<CoordMap>;
  static constexpr Options::String help = {
      "A composition of two UniformTranslations."};

  using options = tmpl::list<
      OptionTags::TimeDependenceCompositionTag<UniformTranslation<MeshDim>>,
      OptionTags::TimeDependenceCompositionTag<UniformTranslation<MeshDim, 1>>>;

  CompositionUniformTranslation() = default;
  ~CompositionUniformTranslation() override = default;
  CompositionUniformTranslation(const CompositionUniformTranslation&) = default;
  CompositionUniformTranslation& operator=(
      const CompositionUniformTranslation&) = default;
  CompositionUniformTranslation(CompositionUniformTranslation&&) = default;
  CompositionUniformTranslation& operator=(CompositionUniformTranslation&&) =
      default;

  explicit CompositionUniformTranslation(
      const std::unique_ptr<TimeDependence<MeshDim>>& uniform_translation0,
      const std::unique_ptr<TimeDependence<MeshDim>>& uniform_translation1);

  auto get_clone() const -> std::unique_ptr<TimeDependence<MeshDim>> override;

  auto block_maps(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, MeshDim>>> override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

 private:
  std::unique_ptr<TimeDependence<MeshDim>> uniform_translation0_;
  std::unique_ptr<TimeDependence<MeshDim>> uniform_translation1_;

  CoordMap coord_map_;
};
}  // namespace time_dependence
}  // namespace creators
}  // namespace domain
