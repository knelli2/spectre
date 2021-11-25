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
#include "Domain/Creators/TimeDependence/GenerateCoordinateMap.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Options/Options.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace domain {
namespace creators {
namespace time_dependence {
/*!
 * \brief A tag used by the `Composition` class to generate a TimeDependence
 * that is a composition of existing `TimeDependences`.
 *
 * The first template parameter is the existing TimeDependence.
 */
template <typename TimeDep>
struct TimeDependenceCompositionTag {
  static constexpr size_t mesh_dim = TimeDep::mesh_dim;
  static std::string name() { return TimeDep::name(); }
  using type = TimeDep;
  static constexpr Options::String help = {
      "One of the maps in the composition."};
  using time_dependence = TimeDep;
};

/*!
 * \brief A TimeDependence that is a composition of various other
 * TimeDependences.
 *
 * To create a new Composition TimeDependence you must create an explicit
 * instantiation of the `Composition<Tags...>` TimeDependence in a `.cpp` file.
 * You must add the `Composition` to the `creatable_classes` list of
 * `TimeDependence` in order for the new Composition to be factory creatable.
 *
 * The tags in the template parameters must be `TimeDependenceCompositionTag`s.
 * See the documentation of the `TimeDependenceCompositionTag` class for details
 * on the tags.
 */
template <typename TimeDependenceCompTag0, typename... TimeDependenceCompTags>
class Composition final
    : public TimeDependence<TimeDependenceCompTag0::mesh_dim> {
 public:
  using CoordMap = detail::generate_coordinate_map_t<tmpl::flatten<
      tmpl::list<typename TimeDependenceCompTag0::time_dependence::
                     MapForComposition::maps_list,
                 typename TimeDependenceCompTags::time_dependence::
                     MapForComposition::maps_list...>>>;
  using maps_list = tmpl::list<CoordMap>;
  static constexpr Options::String help = {"A composition of TimeDependences."};

  static constexpr size_t mesh_dim = TimeDependenceCompTag0::mesh_dim;

  static_assert(
      tmpl::all<
          tmpl::integral_list<size_t, TimeDependenceCompTags::mesh_dim...>,
          std::is_same<tmpl::integral_constant<size_t, mesh_dim>,
                       tmpl::_1>>::value,
      "All TimeDependences passed to Composition must be of the same "
      "dimensionality.");

  using options = tmpl::list<TimeDependenceCompTag0, TimeDependenceCompTags...>;

  Composition() = default;
  ~Composition() override = default;
  Composition(const Composition&) = default;
  Composition& operator=(const Composition&) = default;
  Composition(Composition&&) = default;
  Composition& operator=(Composition&&) = default;

  explicit Composition(
      tmpl::type_from<TimeDependenceCompTag0> first_time_dep,
      tmpl::type_from<TimeDependenceCompTags>... rest_time_dep);

  // static_assert above guarantees that all the TimeDependences have the same
  // dimension so this is ok
  /// Constructor for copying the composition time dependence.
  Composition(
      CoordMap coord_map,
      const std::vector<
          std::unique_ptr<TimeDependence<TimeDependenceCompTag0::mesh_dim>>>&
          time_deps);

  auto get_clone() const -> std::unique_ptr<TimeDependence<mesh_dim>> override;

  auto block_maps(size_t number_of_blocks) const
      -> std::vector<std::unique_ptr<domain::CoordinateMapBase<
          Frame::Grid, Frame::Inertial, mesh_dim>>> override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

 private:
  CoordMap coord_map_;

  // static_assert above guarantees that all the TimeDependences have the same
  // dimension so this is ok
  std::vector<std::unique_ptr<TimeDependence<TimeDependenceCompTag0::mesh_dim>>>
      time_deps_{};
};
}  // namespace time_dependence
}  // namespace creators
}  // namespace domain
