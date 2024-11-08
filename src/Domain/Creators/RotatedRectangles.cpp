// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/RotatedRectangles.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "DataStructures/Index.hpp"
#include "Domain/BoundaryConditions/None.hpp"
#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Domain/Creators/ExpandOverBlocks.hpp"
#include "Domain/Domain.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace domain::creators {
RotatedRectangles::RotatedRectangles(
    const typename LowerBound::type lower_xy,
    const typename Midpoint::type midpoint_xy,
    const typename UpperBound::type upper_xy,
    const typename InitialRefinement::type& initial_refinement_level_xy,
    const typename InitialGridPoints::type initial_number_of_grid_points_in_xy,
    const typename IsPeriodicIn::type is_periodic_in,
    const Options::Context& context)
    // clang-tidy: trivially copyable
    : lower_xy_(std::move(lower_xy)),                       // NOLINT
      midpoint_xy_(std::move(midpoint_xy)),                 // NOLINT
      upper_xy_(std::move(upper_xy)),                       // NOLINT
      is_periodic_in_(std::move(is_periodic_in)),           // NOLINT
      initial_number_of_grid_points_in_xy_(                 // NOLINT
          std::move(initial_number_of_grid_points_in_xy)),  // NOLINT
      boundary_condition_(nullptr) {
  block_names_ = std::vector<std::string>{"LowerLeft", "LowerRight",
                                          "UpperLeft", "UpperRight"};
  block_groups_["All"] =
      std::unordered_set<std::string>{block_names_.begin(), block_names_.end()};

  const ExpandOverBlocks<std::array<size_t, 2>> expand_over_blocks{
      block_names_, block_groups_};
  try {
    initial_refinement_level_xy_ =
        std::visit(expand_over_blocks, initial_refinement_level_xy);
    size_t tmp = 0;
    for (size_t i : {2_st, 3_st}) {
      tmp = initial_refinement_level_xy_[i][0];
      initial_refinement_level_xy_[i][0] = initial_refinement_level_xy_[i][1];
      initial_refinement_level_xy_[i][1] = tmp;
    }
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialRefinement': " << error.what());
  }
}

RotatedRectangles::RotatedRectangles(
    const typename LowerBound::type lower_xy,
    const typename Midpoint::type midpoint_xy,
    const typename UpperBound::type upper_xy,
    const typename InitialRefinement::type& initial_refinement_level_xy,
    const typename InitialGridPoints::type initial_number_of_grid_points_in_xy,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        boundary_condition,
    const Options::Context& context)
    : lower_xy_(lower_xy),
      midpoint_xy_(midpoint_xy),
      upper_xy_(upper_xy),
      is_periodic_in_{{false, false}},
      initial_number_of_grid_points_in_xy_(initial_number_of_grid_points_in_xy),
      boundary_condition_(std::move(boundary_condition)) {
  using domain::BoundaryConditions::is_none;
  if (is_none(boundary_condition_)) {
    PARSE_ERROR(
        context,
        "None boundary condition is not supported. If you would like an "
        "outflow-type boundary condition, you must use that.");
  }
  using domain::BoundaryConditions::is_periodic;
  if (is_periodic(boundary_condition_)) {
    is_periodic_in_[0] = true;
    is_periodic_in_[1] = true;
  }

  block_names_ = std::vector<std::string>{"LowerLeft", "LowerRight",
                                          "UpperLeft", "UpperRight"};
  block_groups_["All"] =
      std::unordered_set<std::string>{block_names_.begin(), block_names_.end()};

  const ExpandOverBlocks<std::array<size_t, 2>> expand_over_blocks{
      block_names_, block_groups_};
  try {
    initial_refinement_level_xy_ =
        std::visit(expand_over_blocks, initial_refinement_level_xy);
    size_t tmp = 0;
    for (size_t i : {2_st, 3_st}) {
      tmp = initial_refinement_level_xy_[i][0];
      initial_refinement_level_xy_[i][0] = initial_refinement_level_xy_[i][1];
      initial_refinement_level_xy_[i][1] = tmp;
    }
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialRefinement': " << error.what());
  }
}

Domain<2> RotatedRectangles::create_domain() const {
  return rectilinear_domain<2>(
      Index<2>{2, 2},
      std::array<std::vector<double>, 2>{
          {{lower_xy_[0], midpoint_xy_[0], upper_xy_[0]},
           {lower_xy_[1], midpoint_xy_[1], upper_xy_[1]}}},
      {},
      std::vector<OrientationMap<2>>{
          OrientationMap<2>::create_aligned(),
          OrientationMap<2>{
              std::array{Direction<2>::lower_xi(), Direction<2>::lower_eta()}},
          OrientationMap<2>{
              std::array{Direction<2>::lower_eta(), Direction<2>::upper_xi()}},
          OrientationMap<2>{
              std::array{Direction<2>::upper_eta(), Direction<2>::lower_xi()}}},
      is_periodic_in_);
}

std::vector<DirectionMap<
    2, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
RotatedRectangles::external_boundary_conditions() const {
  if (boundary_condition_ == nullptr) {
    return {};
  }
  std::vector<DirectionMap<
      2, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
      boundary_conditions{4};
  if (is_periodic_in_[0]) {
    return boundary_conditions;
  }
  boundary_conditions[0][Direction<2>::lower_xi()] =
      boundary_condition_->get_clone();
  boundary_conditions[0][Direction<2>::lower_eta()] =
      boundary_condition_->get_clone();
  boundary_conditions[1][Direction<2>::lower_xi()] =
      boundary_condition_->get_clone();
  boundary_conditions[1][Direction<2>::upper_eta()] =
      boundary_condition_->get_clone();
  boundary_conditions[2][Direction<2>::upper_xi()] =
      boundary_condition_->get_clone();
  boundary_conditions[2][Direction<2>::upper_eta()] =
      boundary_condition_->get_clone();
  boundary_conditions[3][Direction<2>::lower_xi()] =
      boundary_condition_->get_clone();
  boundary_conditions[3][Direction<2>::upper_eta()] =
      boundary_condition_->get_clone();
  return boundary_conditions;
}

std::vector<std::array<size_t, 2>> RotatedRectangles::initial_extents() const {
  const size_t& x_0 = initial_number_of_grid_points_in_xy_[0][0];
  const size_t& x_1 = initial_number_of_grid_points_in_xy_[0][1];
  const size_t& y_0 = initial_number_of_grid_points_in_xy_[1][0];
  const size_t& y_1 = initial_number_of_grid_points_in_xy_[1][1];
  return {{{x_0, y_0}}, {{x_1, y_0}}, {{y_1, x_0}}, {{y_1, x_1}}};
}

std::vector<std::array<size_t, 2>>
RotatedRectangles::initial_refinement_levels() const {
  return initial_refinement_level_xy_;
}
}  // namespace domain::creators
