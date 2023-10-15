// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class RotatedIntervals.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

#include "Domain/BoundaryConditions/GetBoundaryConditionsBase.hpp"
#include "Domain/Creators/DomainCreator.hpp"  // IWYU pragma: keep
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain {
namespace CoordinateMaps {
class Affine;
template <size_t VolumeDim>
class DiscreteRotation;
}  // namespace CoordinateMaps

template <typename SourceFrame, typename TargetFrame, typename... Maps>
class CoordinateMap;
}  // namespace domain
/// \endcond

namespace domain {
namespace creators {
/// Create a 1D Domain consisting of two rotated Blocks.
/// The left block has its logical \f$\xi\f$-axis aligned with the grid x-axis.
/// The right block has its logical \f$\xi\f$-axis opposite to the grid x-axis.
/// This is useful for testing code that deals with unaligned blocks.
class RotatedIntervals : public DomainCreator<1> {
 public:
  using maps_list = tmpl::list<domain::CoordinateMap<
      Frame::BlockLogical, Frame::Inertial,
      domain::CoordinateMaps::DiscreteRotation<1>, CoordinateMaps::Affine>>;

  struct LowerBound {
    using type = std::array<double, 1>;
    inline const static std::string help {
        "Sequence of [x], the lower bound in the target frame."};
  };

  struct Midpoint {
    using type = std::array<double, 1>;
    inline const static std::string help {
        "Sequence of [x], the midpoint in the target frame."};
  };

  struct UpperBound {
    using type = std::array<double, 1>;
    inline const static std::string help {
        "Sequence of [x], the upper bound in the target frame."};
  };

  struct IsPeriodicIn {
    using type = std::array<bool, 1>;
    inline const static std::string help {
        "Sequence for [x], true if periodic."};
  };
  struct InitialRefinement {
    using type = std::array<size_t, 1>;
    inline const static std::string help {
        "Initial refinement level in [x]."};
  };

  struct InitialGridPoints {
    using type = std::array<std::array<size_t, 2>, 1>;
    inline const static std::string help {
        "Initial number of grid points in [[x]]."};
  };

  struct TimeDependence {
    using type =
        std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>;
    inline const static std::string help {
        "The time dependence of the moving mesh domain."};
  };

  struct BoundaryConditions {
    inline const static std::string help {"The boundary conditions to apply."};
  };
  template <typename BoundaryConditionsBase>
  struct UpperBoundaryCondition {
    static std::string name() { return "UpperBoundary"; }
    inline const static std::string help
        {"Options for the boundary condition applied at the upper boundary."};
    using type = std::unique_ptr<BoundaryConditionsBase>;
    using group = BoundaryConditions;
  };
  template <typename BoundaryConditionsBase>
  struct LowerBoundaryCondition {
    static std::string name() { return "LowerBoundary"; }
    inline const static std::string help
        {"Options for the boundary condition applied at the lower boundary."};
    using type = std::unique_ptr<BoundaryConditionsBase>;
    using group = BoundaryConditions;
  };

  using common_options = tmpl::list<LowerBound, Midpoint, UpperBound,
                                    InitialRefinement, InitialGridPoints>;
  using options_periodic = tmpl::list<IsPeriodicIn>;
  template <typename System>
  using options_boundary_conditions = tmpl::list<
      LowerBoundaryCondition<
          domain::BoundaryConditions::get_boundary_conditions_base<System>>,
      UpperBoundaryCondition<
          domain::BoundaryConditions::get_boundary_conditions_base<System>>>;

  template <typename Metavariables>
  using options = tmpl::append<
      common_options,
      tmpl::conditional_t<
          domain::BoundaryConditions::has_boundary_conditions_base_v<
              typename Metavariables::system>,
          options_boundary_conditions<typename Metavariables::system>,
          options_periodic>,
      tmpl::list<TimeDependence>>;

  inline const static std::string help {
      "A DomainCreator useful for testing purposes.\n"
      "RotatedIntervals creates the interval [LowerX,UpperX] from two\n"
      "rotated Blocks. The outermost index to InitialGridPoints is the\n"
      "dimension index (of which there is only one in the case of\n"
      "RotatedIntervals), and the innermost index is the block index\n"
      "along that dimension."};

  RotatedIntervals(
      std::array<double, 1> lower_x, std::array<double, 1> midpoint_x,
      std::array<double, 1> upper_x,
      std::array<size_t, 1> initial_refinement_level_x,
      std::array<std::array<size_t, 2>, 1> initial_number_of_grid_points_in_x,
      std::array<bool, 1> is_periodic_in,
      std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
          time_dependence);

  RotatedIntervals(
      std::array<double, 1> lower_x, std::array<double, 1> midpoint_x,
      std::array<double, 1> upper_x,
      std::array<size_t, 1> initial_refinement_level_x,
      std::array<std::array<size_t, 2>, 1> initial_number_of_grid_points_in_x,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          lower_boundary_condition,
      std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
          upper_boundary_condition,
      std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
          time_dependence,
      const Options::Context& context = {});

  RotatedIntervals() = default;
  RotatedIntervals(const RotatedIntervals&) = delete;
  RotatedIntervals(RotatedIntervals&&) = default;
  RotatedIntervals& operator=(const RotatedIntervals&) = delete;
  RotatedIntervals& operator=(RotatedIntervals&&) = default;
  ~RotatedIntervals() override = default;

  Domain<1> create_domain() const override;

  std::vector<DirectionMap<
      1, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
  external_boundary_conditions() const override;

  std::vector<std::array<size_t, 1>> initial_extents() const override;

  std::vector<std::array<size_t, 1>> initial_refinement_levels() const override;

  auto functions_of_time(const std::unordered_map<std::string, double>&
                             initial_expiration_times = {}) const
      -> std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>> override;

 private:
  std::array<double, 1> lower_x_{
      {std::numeric_limits<double>::signaling_NaN()}};
  std::array<double, 1> midpoint_x_{
      {std::numeric_limits<double>::signaling_NaN()}};
  std::array<double, 1> upper_x_{
      {std::numeric_limits<double>::signaling_NaN()}};
  std::array<bool, 1> is_periodic_in_{{false}};
  std::array<size_t, 1> initial_refinement_level_x_{
      {std::numeric_limits<size_t>::max()}};
  std::array<std::array<size_t, 2>, 1> initial_number_of_grid_points_in_x_{
      {{{std::numeric_limits<size_t>::max()}}}};
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      lower_boundary_condition_{nullptr};
  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
      upper_boundary_condition_{nullptr};
  std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
      time_dependence_{nullptr};
};
}  // namespace creators
}  // namespace domain
