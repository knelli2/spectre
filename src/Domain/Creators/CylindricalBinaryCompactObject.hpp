// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/BoundaryConditions/GetBoundaryConditionsBase.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/Creators/BinaryCompactObjectHelpers.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace CoordinateMaps {
class Interval;
template <typename Map1, typename Map2>
class ProductOf2Maps;
template <typename Map1, typename Map2, typename Map3>
class ProductOf3Maps;
template <size_t VolumeDim>
class Wedge;
template <size_t VolumeDim>
class DiscreteRotation;
class UniformCylindricalEndcap;
class UniformCylindricalFlatEndcap;
class UniformCylindricalSide;
}  // namespace CoordinateMaps

template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;

namespace FunctionsOfTime {
class FunctionOfTime;
}  // namespace FunctionsOfTime
}  // namespace domain

namespace Frame {
struct Grid;
struct Distorted;
struct Inertial;
struct BlockLogical;
}  // namespace Frame
/// \endcond

namespace domain::creators {

/*!
 * \ingroup ComputationalDomainGroup
 *
 * \brief A general domain for two compact objects based on cylinders.
 *
 * Creates a 3D Domain that represents a binary compact object
 * solution.  This domain is described briefly in the Appendix of
 * \cite Buchman:2012dw, and is illustrated in Figure 20 of that
 * paper.
 *
 * In the code and options below, `ObjectA` and `ObjectB` refer to the
 * two compact objects. In the grid frame, `ObjectA` is located to the
 * right of (i.e. a more positive value of the x-coordinate than)
 * `ObjectB`.  The inner edge of the Blocks surrounding each of
 * `ObjectA` and `ObjectB` is spherical in grid coordinates; the
 * user must specify the center and radius of this surface for both
 * `ObjectA` and `ObjectB`, and the user must specify the outer boundary
 * radius.  The outer boundary is a sphere centered at the origin.
 *
 * This domain offers some grid anchors. See
 * `domain::creators::bco::create_grid_anchors` for which ones are offered.
 *
 * Note that Figure 20 of \cite Buchman:2012dw illustrates additional
 * spherical shells inside the "EA" and "EB" blocks, and the caption
 * of Figure 20 indicates that there are additional spherical shells
 * outside the "CA" and "CB" blocks; `CylindricalBinaryCompactObject`
 * has these extra shells inside "EA" only if the option `IncludeInnerSphereA`
 * is true, it has the extra shells inside "EB" only if the option
 * `IncludeInnerSphereB` is true, and it has the extra shells outside
 * "CA" and "CB" only if `IncludeOuterSphere` is true.
 * If the shells are absent, then the "EA" and "EB"
 * blocks extend to the excision boundaries and the "CA" and "CB" blocks
 * extend to the outer boundary.
 *
 * The Blocks are named as follows:
 * - Each of CAFilledCylinder, EAFilledCylinder, EBFilledCylinder,
 *   MAFilledCylinder, MBFilledCylinder, and CBFilledCylinder consists
 *   of 5 blocks, named 'Center', 'East', 'North', 'West', and
 *   'South', so an example of a valid block name is
 *   'CAFilledCylinderCenter'.
 * - Each of CACylinder, EACylinder, EBCylinder, and CBCylinder
 *   consists of 4 blocks, named 'East', 'North', 'West', and 'South',
 *   so an example of a valid block name is 'CACylinderEast'.
 * - The Block group called "Outer" consists of all the CA and CB blocks. They
 *   all border the outer boundary if `IncludeOuterSphere` is false.
 * - If `IncludeOuterSphere` is true, then there are more blocks named
 *   OuterSphereCAFilledCylinder, OuterSphereCBFilledCylinder,
 *   OuterSphereCACylinder, and OuterSphereCBCylinder.
 *   These are in a Block group called "OuterSphere",
 *   and all of these border the outer boundary.
 * - The Block group called "InnerA" consists of all the EA, and MA
 *   blocks. They all border the inner boundary "A" if
 *   `IncludeInnerSphereA` is false.
 * - If `IncludeInnerSphereA` is true, then there are new blocks
 *   InnerSphereEAFilledCylinder, InnerSphereMAFilledCylinder, and
 *   InnerSphereEACylinder. These are in a Block group called "InnerSphereA",
 *   and all of these border the inner excision boundary "A".
 * - The Block group called "InnerB" consists of all the EB, and MB
 *   blocks. They all border the inner boundary "B" if
 *   `IncludeInnerSphereB` is false.
 * - If `IncludeInnerSphereB` is true, then there are new blocks
 *   InnerSphereEBFilledCylinder, InnerSphereMBFilledCylinder, and
 *   InnerSphereEBCylinder. These are in a Block group called "InnerSphereB",
 *   and all of these border the inner excision boundary "B".
 *
 * If \f$c_A\f$ and \f$c_B\f$ are the input parameters center_A and
 * center_B, \f$r_A\f$ and \f$r_B\f$ are the input parameters radius_A and
 * radius_B, and \f$R\f$ is the outer boundary radius, we demand the
 * following restrictions on parameters:
 * - \f$c_A^0>0\f$; this is a convention to simplify the code.
 * - \f$c_B^0<0\f$; this is a convention to simplify the code.
 * - \f$|c_A^0|\le|c_B^0|\f$. We should roughly have \f$r_A c_A^0 + r_B c_B^0\f$
 *   close to zero; that is, for BBHs (where \f$r_A\f$ is roughly twice the
 *   mass of the heavier object A, and \f$r_B\f$ is roughly twice the mass
 *   of the lighter object B) the center of mass should be roughly
 *   at the origin.
 * - \f$0 < r_B < r_A\f$
 * - \f$R \ge 3(|c_A^0|-|c_B^0|)\f$; otherwise the blocks will be too compressed
 *   near the outer boundary.
 *
 * All time dependent maps are optional to specify. To include a map, specify
 * its options. Otherwise specify `None` for that map. You can also turn off
 * time dependent maps all together by specifying `None` for the
 * `TimeDependentMaps` option. See
 * `domain::creators::bco::TimeDependentMapOptions`.
 */
class CylindricalBinaryCompactObject : public DomainCreator<3> {
 public:
  using maps_list = tmpl::flatten<
      tmpl::list<domain::CoordinateMap<
                     Frame::BlockLogical, Frame::Inertial,
                     CoordinateMaps::ProductOf3Maps<CoordinateMaps::Interval,
                                                    CoordinateMaps::Interval,
                                                    CoordinateMaps::Interval>,
                     CoordinateMaps::UniformCylindricalEndcap,
                     CoordinateMaps::DiscreteRotation<3>>,
                 domain::CoordinateMap<
                     Frame::BlockLogical, Frame::Inertial,
                     CoordinateMaps::ProductOf2Maps<CoordinateMaps::Wedge<2>,
                                                    CoordinateMaps::Interval>,
                     CoordinateMaps::UniformCylindricalEndcap,
                     CoordinateMaps::DiscreteRotation<3>>,
                 domain::CoordinateMap<
                     Frame::BlockLogical, Frame::Inertial,
                     CoordinateMaps::ProductOf3Maps<CoordinateMaps::Interval,
                                                    CoordinateMaps::Interval,
                                                    CoordinateMaps::Interval>,
                     CoordinateMaps::UniformCylindricalFlatEndcap,
                     CoordinateMaps::DiscreteRotation<3>>,
                 domain::CoordinateMap<
                     Frame::BlockLogical, Frame::Inertial,
                     CoordinateMaps::ProductOf2Maps<CoordinateMaps::Wedge<2>,
                                                    CoordinateMaps::Interval>,
                     CoordinateMaps::UniformCylindricalFlatEndcap,
                     CoordinateMaps::DiscreteRotation<3>>,
                 domain::CoordinateMap<
                     Frame::BlockLogical, Frame::Inertial,
                     CoordinateMaps::ProductOf2Maps<CoordinateMaps::Wedge<2>,
                                                    CoordinateMaps::Interval>,
                     CoordinateMaps::UniformCylindricalSide,
                     CoordinateMaps::DiscreteRotation<3>>,
                 bco::TimeDependentMapOptions::maps_list>>;

  struct CenterA {
    using type = std::array<double, 3>;
    inline const static std::string help {
        "Grid coordinates of center for Object A, which is at x>0."};
  };
  struct CenterB {
    using type = std::array<double, 3>;
    inline const static std::string help {
        "Grid coordinates of center for Object B, which is at x<0."};
  };
  struct RadiusA {
    using type = double;
    inline const static std::string help {
        "Grid-coordinate radius of grid boundary around Object A."};
  };
  struct RadiusB {
    using type = double;
    inline const static std::string help {
        "Grid-coordinate radius of grid boundary around Object B."};
  };
  struct IncludeInnerSphereA {
    using type = bool;
    inline const static std::string help {
        "Add an extra spherical layer of Blocks around Object A."};
  };
  struct IncludeInnerSphereB {
    using type = bool;
    inline const static std::string help {
        "Add an extra spherical layer of Blocks around Object B."};
  };
  struct IncludeOuterSphere {
    using type = bool;
    inline const static std::string help {
        "Add an extra spherical layer of Blocks inside the outer boundary."};
  };
  struct OuterRadius {
    using type = double;
    inline const static std::string help {
        "Grid-coordinate radius of outer boundary."};
  };
  struct UseEquiangularMap {
    using type = bool;
    inline const static std::string help {
        "Distribute grid points equiangularly in 2d wedges."};
    static bool suggested_value() { return false; }
  };

  struct InitialRefinement {
    using type =
        std::variant<size_t, std::array<size_t, 3>,
                     std::vector<std::array<size_t, 3>>,
                     std::unordered_map<std::string, std::array<size_t, 3>>>;
    inline const static std::string help {
        "Initial refinement level. Specify one of: a single number, a list "
        "representing [r, theta, perp], or such a list for every block in the "
        "domain. Here 'r' is the radial direction normal to the inner and "
        "outer boundaries, 'theta' is the periodic direction, and 'perp' is "
        "the third direction."};
  };
  struct InitialGridPoints {
    using type =
        std::variant<size_t, std::array<size_t, 3>,
                     std::vector<std::array<size_t, 3>>,
                     std::unordered_map<std::string, std::array<size_t, 3>>>;
    inline const static std::string help {
        "Initial number of grid points. Specify one of: a single number, a "
        "list representing [r, theta, perp], or such a list for every block in "
        "the domain. Here 'r' is the radial direction normal to the inner and "
        "outer boundaries, 'theta' is the periodic direction, and 'perp' is "
        "the third direction."};
  };

  struct BoundaryConditions {
    inline const static std::string help {"The boundary conditions to apply."};
  };
  template <typename BoundaryConditionsBase>
  struct InnerBoundaryCondition {
    static std::string name() { return "InnerBoundary"; }
    inline const static std::string help
        {"Options for the inner boundary conditions."};
    using type = std::unique_ptr<BoundaryConditionsBase>;
    using group = BoundaryConditions;
  };

  template <typename BoundaryConditionsBase>
  struct OuterBoundaryCondition {
    static std::string name() { return "OuterBoundary"; }
    inline const static std::string help
        {"Options for the outer boundary conditions."};
    using type = std::unique_ptr<BoundaryConditionsBase>;
    using group = BoundaryConditions;
  };

  struct TimeDependentMaps {
    using type =
        Options::Auto<bco::TimeDependentMapOptions, Options::AutoLabel::None>;
    inline const static std::string help =bco::TimeDependentMapOptions::help;
  };

  using time_independent_options =
      tmpl::list<CenterA, CenterB, RadiusA, RadiusB, IncludeInnerSphereA,
                 IncludeInnerSphereB, IncludeOuterSphere, OuterRadius,
                 UseEquiangularMap, InitialRefinement, InitialGridPoints>;

  template <typename Metavariables>
  using basic_options = tmpl::conditional_t<
      domain::creators::bco::enable_time_dependent_maps_v<Metavariables>,
      tmpl::push_front<time_independent_options, TimeDependentMaps>,
      time_independent_options>;

  template <typename Metavariables>
  using options = tmpl::conditional_t<
      domain::BoundaryConditions::has_boundary_conditions_base_v<
          typename Metavariables::system>,
      tmpl::push_back<
          basic_options<Metavariables>,
          InnerBoundaryCondition<
              domain::BoundaryConditions::get_boundary_conditions_base<
                  typename Metavariables::system>>,
          OuterBoundaryCondition<
              domain::BoundaryConditions::get_boundary_conditions_base<
                  typename Metavariables::system>>>,
      basic_options<Metavariables>>;

  inline const static std::string help{
      "The CylindricalBinaryCompactObject domain is a general domain for "
      "two compact objects. The user must provide the (grid-frame) "
      "centers and radii of the spherical inner edge of the grid surrounding "
      "each of the two compact objects A and B."};

  CylindricalBinaryCompactObject(
      std::array<double, 3> center_A, std::array<double, 3> center_B,
      double radius_A, double radius_B, bool include_inner_sphere_A,
      bool include_inner_sphere_B, bool include_outer_sphere,
      double outer_radius, bool use_equiangular_map,
      const typename InitialRefinement::type& initial_refinement,
      const typename InitialGridPoints::type& initial_grid_points,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          inner_boundary_condition = nullptr,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          outer_boundary_condition = nullptr,
      const Options::Context& context = {});

  CylindricalBinaryCompactObject(
      std::optional<bco::TimeDependentMapOptions> time_dependent_options,
      std::array<double, 3> center_A, std::array<double, 3> center_B,
      double radius_A, double radius_B, bool include_inner_sphere_A,
      bool include_inner_sphere_B, bool include_outer_sphere,
      double outer_radius, bool use_equiangular_map,
      const typename InitialRefinement::type& initial_refinement,
      const typename InitialGridPoints::type& initial_grid_points,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          inner_boundary_condition = nullptr,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          outer_boundary_condition = nullptr,
      const Options::Context& context = {});

  CylindricalBinaryCompactObject() = default;
  CylindricalBinaryCompactObject(const CylindricalBinaryCompactObject&) =
      delete;
  CylindricalBinaryCompactObject(CylindricalBinaryCompactObject&&) = default;
  CylindricalBinaryCompactObject& operator=(
      const CylindricalBinaryCompactObject&) = delete;
  CylindricalBinaryCompactObject& operator=(CylindricalBinaryCompactObject&&) =
      default;
  ~CylindricalBinaryCompactObject() override = default;

  Domain<3> create_domain() const override;

  std::unordered_map<std::string, tnsr::I<double, 3, Frame::Grid>>
  grid_anchors() const override {
    return grid_anchors_;
  }

  std::vector<DirectionMap<
      3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
  external_boundary_conditions() const override;

  std::vector<std::array<size_t, 3>> initial_extents() const override;

  std::vector<std::array<size_t, 3>> initial_refinement_levels() const override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

  std::vector<std::string> block_names() const override { return block_names_; }

  std::unordered_map<std::string, std::unordered_set<std::string>>
  block_groups() const override {
    return block_groups_;
  }

 private:
  // Note that center_A_ and center_B_ are rotated with respect to the
  // input centers (which are in the grid frame), so that we can
  // construct the map in a frame where the centers are offset in the
  // z direction.  At the end, there will be another rotation back to
  // the grid frame (where the centers are offset in the x direction).
  std::array<double, 3> center_A_{};
  std::array<double, 3> center_B_{};
  double radius_A_{};
  double radius_B_{};
  double outer_radius_A_{};
  double outer_radius_B_{};
  bool include_inner_sphere_A_{};
  bool include_inner_sphere_B_{};
  bool include_outer_sphere_{};
  double outer_radius_{};
  bool use_equiangular_map_{false};
  typename std::vector<std::array<size_t, 3>> initial_refinement_{};
  typename std::vector<std::array<size_t, 3>> initial_grid_points_{};
  // cut_spheres_offset_factor_ is eta in Eq. (A.9) of
  // https://arxiv.org/abs/1206.3015.  cut_spheres_offset_factor_
  // could be set to unity to simplify the equations.  Here we fix it
  // to the value 0.99 used in SpEC, so that we reproduce SpEC's
  // domain decomposition.
  double cut_spheres_offset_factor_{0.99};
  // z_cutting_plane_ is x_C in Eq. (A.9) of
  // https://arxiv.org/abs/1206.3015 (but rotated to the z-axis).
  double z_cutting_plane_{};
  size_t number_of_blocks_{};
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      inner_boundary_condition_;
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      outer_boundary_condition_;
  std::vector<std::string> block_names_{};
  std::unordered_map<std::string, std::unordered_set<std::string>>
      block_groups_{};
  std::unordered_map<std::string, tnsr::I<double, 3, Frame::Grid>>
      grid_anchors_{};
  // FunctionsOfTime options
  std::optional<bco::TimeDependentMapOptions> time_dependent_options_{};
};
}  // namespace domain::creators
