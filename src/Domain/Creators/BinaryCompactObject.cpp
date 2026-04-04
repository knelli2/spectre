// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/BinaryCompactObject.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Block.hpp"
#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/CoordinateMaps/Equiangular.hpp"
#include "Domain/CoordinateMaps/Frustum.hpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/Interval.hpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/CoordinateMaps/SphericalToCartesianPfaffian.hpp"
#include "Domain/CoordinateMaps/Wedge.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/ExpandOverBlocks.hpp"
#include "Domain/Creators/ShellDistribution.hpp"
#include "Domain/Creators/TimeDependentOptions/BinaryCompactObject.hpp"
#include "Domain/Domain.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/ExcisionSphere.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/Structure/BlockNeighbors.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/Topology.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"

namespace Frame {
struct BlockLogical;
}  // namespace Frame

namespace domain::creators {
namespace bco {
std::unordered_map<std::string, tnsr::I<double, 3, Frame::Grid>>
create_grid_anchors(const std::array<double, 3>& center_a,
                    const std::array<double, 3>& center_b) {
  std::unordered_map<std::string, tnsr::I<double, 3, Frame::Grid>> result{};
  result["Center" + get_output(ObjectLabel::A)] =
      tnsr::I<double, 3, Frame::Grid>{center_a};
  result["Center" + get_output(ObjectLabel::B)] =
      tnsr::I<double, 3, Frame::Grid>{center_b};
  result["Center"] =
      tnsr::I<double, 3, Frame::Grid>{0.5 * (center_a + center_b)};

  return result;
}
}  // namespace bco

template <bool UseWorldtube>
bool BinaryCompactObject<UseWorldtube>::Object::is_excised() const {
  return inner_boundary_condition.has_value();
}

template <bool UseWorldtube>
BinaryCompactObject<UseWorldtube>::BinaryCompactObject(
    typename ObjectA::type object_A, typename ObjectB::type object_B,
    std::array<double, 2> center_of_mass_offset, const double envelope_radius,
    const double outer_radius, const double cube_scale,
    const typename InitialRefinement::type& initial_refinement,
    const typename InitialGridPoints::type& initial_number_of_grid_points,
    const bool use_equiangular_map,
    const CoordinateMaps::Distribution radial_distribution_envelope,
    const std::vector<double>& radial_partitioning_outer_shell,
    const typename RadialDistributionOuterShell::type&
        radial_distribution_outer_shell,
    const double opening_angle_in_degrees,
    std::optional<bco::TimeDependentMapOptions<false>> time_dependent_options,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        outer_boundary_condition,
    const Options::Context& context)
    : object_A_(std::move(object_A)),
      object_B_(std::move(object_B)),
      center_of_mass_offset_(center_of_mass_offset),
      envelope_radius_(envelope_radius),
      outer_radius_(outer_radius),
      use_equiangular_map_(use_equiangular_map),
      radial_distribution_envelope_(radial_distribution_envelope),
      radial_partitioning_outer_shell_(radial_partitioning_outer_shell),
      outer_boundary_condition_(std::move(outer_boundary_condition)),
      x_coord_a_(
          std::visit([](const auto& arg) { return arg.x_coord; }, object_A_)),
      x_coord_b_(
          std::visit([](const auto& arg) { return arg.x_coord; }, object_B_)),
      is_excised_a_(std::visit([](const auto& arg) { return arg.is_excised(); },
                               object_A_)),
      is_excised_b_(std::visit([](const auto& arg) { return arg.is_excised(); },
                               object_B_)),
      use_single_block_a_(
          std::holds_alternative<CartesianCubeAtXCoord>(object_A_)),
      use_single_block_b_(
          std::holds_alternative<CartesianCubeAtXCoord>(object_B_)),
      time_dependent_options_(std::move(time_dependent_options)),
      opening_angle_(M_PI * opening_angle_in_degrees / 180.0) {
  // Determination of parameters for domain construction:
  const double tan_half_opening_angle = tan(0.5 * opening_angle_);
  translation_ = 0.5 * (x_coord_a_ + x_coord_b_);
  length_inner_cube_ = cube_scale * (x_coord_a_ - x_coord_b_);
  if (length_inner_cube_ < (x_coord_a_ - x_coord_b_)) {
    PARSE_ERROR(
        context,
        "The cube length should be greater than or equal to the initial "
        "separation between the two objects.");
  }
  length_outer_cube_ =
      2.0 * envelope_radius_ / sqrt(2.0 + square(tan_half_opening_angle));

  // We chose to handle potential roundoff differences here by using equal
  // within roundoff instead of exact equality because the wedge map expects
  // exact equality when checking for a zero offset.
  if (equal_within_roundoff(cube_scale, 1.0)) {
    offset_x_coord_a_ = 0.0;
    offset_x_coord_b_ = 0.0;
  } else {
    offset_x_coord_a_ =
        x_coord_a_ - (x_coord_a_ + x_coord_b_ + length_inner_cube_) * 0.5;
    offset_x_coord_b_ =
        x_coord_b_ - (x_coord_a_ + x_coord_b_ - length_inner_cube_) * 0.5;
  }

  set_shell_distribution(make_not_null(&number_of_outer_shells_),
                         make_not_null(&radial_distribution_outer_shell_),
                         radial_partitioning_outer_shell_,
                         radial_distribution_outer_shell, envelope_radius_,
                         outer_radius_, "envelope", "outer", context);

  // Calculate number of blocks
  // Object A's shell has 6 wedge blocks and cube has 6 wedge blocks (12 total).
  // Object B's shell is now a single spherical shell block and cube has 6 wedge
  // blocks (7 total). The envelope and each outer shell have another 10 blocks.
  // Total base: 12 + 7 + 10 = 29.
  number_of_blocks_ = 29 + 10 * number_of_outer_shells_;
  // For each object whose interior is not excised, add 1 block
  if ((not use_single_block_a_) and (not is_excised_a_)) {
    number_of_blocks_++;
  }
  if ((not use_single_block_b_) and (not is_excised_b_)) {
    number_of_blocks_++;
  }

  // For object A replaced by a single block, remove (12-1)=11
  if (use_single_block_a_) {
    number_of_blocks_ -= 11;
  }
  // For object B replaced by a single block, remove (7-1)=6 since B's shell is
  // now 1 spherical shell block and its cube is 6 wedge blocks (7 total).
  if (use_single_block_b_) {
    number_of_blocks_ -= 6;
  }

  if (x_coord_a_ <= 0.0) {
    PARSE_ERROR(
        context,
        "The x-coordinate of ObjectA's center is expected to be positive.");
  }
  if (x_coord_b_ >= 0.0) {
    PARSE_ERROR(
        context,
        "The x-coordinate of ObjectB's center is expected to be negative.");
  }
  if (envelope_radius_ <= length_inner_cube_ * sqrt(3.0)) {
    const double suggested_value = 2.0 * length_inner_cube_ * sqrt(3.0);
    PARSE_ERROR(
        context,
        "The radius for the enveloping cube is too small! The Frustums will be "
        "malformed. A recommended radius is:\n"
            << suggested_value);
  }
  // The following options are irrelevant if the inner regions are covered
  // with simple blocks, so we only check them if object_A_ uses the first
  // variant type.
  if (not use_single_block_a_) {
    const auto& object_a = std::get<Object>(object_A_);
    if (object_a.outer_radius < object_a.inner_radius) {
      PARSE_ERROR(context,
                  "ObjectA's inner radius must be less than its outer radius.");
    }
    // The length of the cube surrounding each object is length_inner_cube_
    // / 2.0, this checks to make sure that the shell is not touching the cube.
    if (object_a.outer_radius >=
        (length_inner_cube_ / 2.0) + offset_x_coord_a_) {
      PARSE_ERROR(
          context,
          "ObjectA's outer radius is too large for the given separation, "
          " try using "
              << length_inner_cube_ / (2.5 * cube_scale));
    }
    if (object_a.use_logarithmic_map and not object_a.is_excised()) {
      PARSE_ERROR(
          context,
          "Using a logarithmically spaced radial grid in the part "
          "of Layer 1 enveloping Object A requires excising the interior of "
          "Object A");
    }
    if (not object_a.is_excised() and
        not equal_within_roundoff(cube_scale, 1.0)) {
      ERROR(
          "A filled object cannot be offset within its cube. Set the CubeScale "
          "to 1.");
    }
    if (object_a.is_excised() and
        ((*object_a.inner_boundary_condition == nullptr) !=
         (outer_boundary_condition_ == nullptr))) {
      PARSE_ERROR(
          context,
          "Must specify either both inner and outer boundary conditions "
          "or neither.");
    }
  }
  if (not use_single_block_b_) {
    const auto& object_b = std::get<Object>(object_B_);
    if (object_b.outer_radius < object_b.inner_radius) {
      PARSE_ERROR(context,
                  "ObjectB's inner radius must be less than its outer radius.");
    }
    // The length of the cube surrounding each object is length_inner_cube_
    // / 2.0, this checks to make sure that the shell is not touching the cube.
    if (object_b.outer_radius >=
        (length_inner_cube_ / 2.0) - offset_x_coord_b_) {
      PARSE_ERROR(
          context,
          "ObjectB's outer radius is too large for the given separation, "
          " try using "
              << length_inner_cube_ / (2.5 * cube_scale));
    }
    if (object_b.use_logarithmic_map and not object_b.is_excised()) {
      PARSE_ERROR(
          context,
          "Using a logarithmically spaced radial grid in the part "
          "of Layer 1 enveloping Object B requires excising the interior of "
          "Object B");
    }
    if (not object_b.is_excised() and
        not equal_within_roundoff(cube_scale, 1.0)) {
      ERROR(
          "A filled object cannot be offset within its cube. Set the CubeScale "
          "to 1.");
    }
    if (object_b.is_excised() and
        ((*object_b.inner_boundary_condition == nullptr) !=
         (outer_boundary_condition_ == nullptr))) {
      PARSE_ERROR(
          context,
          "Must specify either both inner and outer boundary conditions "
          "or neither.");
    }
  }
  const bool filled_excision_a = not(use_single_block_a_ or is_excised_a_);
  const bool filled_excision_b = not(use_single_block_b_ or is_excised_b_);
  if ((filled_excision_a or filled_excision_b) and
      not equal_within_roundoff(offset_x_coord_a_, 0.0)) {
    PARSE_ERROR(context,
                "Setting CubeScale > 1.0 is not supported for domains with "
                "ExciseInterior = False. Consider using "
                "CartesianCubeAtX for the Object without an excised interior.");
  }

  if (envelope_radius_ >= outer_radius_) {
    PARSE_ERROR(context,
                "The outer radius must be larger than the envelope radius.");
  }
  using domain::BoundaryConditions::is_periodic;
  if (is_periodic(outer_boundary_condition_) or
      (is_excised_a_ and
       is_periodic(*(std::get<Object>(object_A_).inner_boundary_condition))) or
      (is_excised_b_ and
       is_periodic(*(std::get<Object>(object_B_).inner_boundary_condition)))) {
    PARSE_ERROR(
        context,
        "Cannot have periodic boundary conditions with a binary domain");
  }

  // Create grid anchors
  grid_anchors_ =
      bco::create_grid_anchors(std::array{x_coord_a_, center_of_mass_offset_[0],
                                          center_of_mass_offset_[1]},
                               std::array{x_coord_b_, center_of_mass_offset_[0],
                                          center_of_mass_offset_[1]});

  // Create block names and groups
  static std::array<std::string, 6> wedge_directions{
      "UpperZ", "LowerZ", "UpperY", "LowerY", "UpperX", "LowerX"};
  const auto add_object_region = [this](const std::string& object_name,
                                        const std::string& region_name) {
    for (const std::string& wedge_direction : wedge_directions) {
      const std::string block_name =
          // NOLINTNEXTLINE(performance-inefficient-string-concatenation)
          object_name + region_name + wedge_direction;
      block_names_.push_back(block_name);
      block_groups_[object_name + region_name].insert(block_name);
    }
  };
  const auto add_object_interior = [this](const std::string& object_name) {
    const std::string block_name = object_name + "Interior";
    block_names_.push_back(block_name);
  };
  const auto add_outer_region = [this](const std::string& region_name) {
    for (const std::string& wedge_direction : wedge_directions) {
      for (const std::string& leftright : {"Left"s, "Right"s}) {
        if ((wedge_direction == "UpperX" and leftright == "Left") or
            (wedge_direction == "LowerX" and leftright == "Right")) {
          // The outer regions are divided in half perpendicular to the
          // x-axis at x=0. Therefore, the left side only has a block in
          // negative x-direction, and the right side only has one in
          // positive x-direction.
          continue;
        }
        // NOLINTNEXTLINE(performance-inefficient-string-concatenation)
        const std::string block_name =
            region_name + wedge_direction +
            (wedge_direction == "UpperX" or wedge_direction == "LowerX"
                 ? ""
                 : leftright);
        block_names_.push_back(block_name);
        block_groups_[region_name].insert(block_name);
      }
    }
  };
  // Finding the first block of outer shell
  first_outer_shell_block_ = 0;
  if (use_single_block_a_) {
    block_names_.emplace_back("ObjectA");
    first_outer_shell_block_ += 1;
  } else {
    add_object_region("ObjectA", "Shell");  // 6 blocks
    add_object_region("ObjectA", "Cube");   // 6 blocks
    first_outer_shell_block_ += 12;
  }
  if (use_single_block_b_) {
    block_names_.emplace_back("ObjectB");
    first_outer_shell_block_ += 1;
  } else {
    // Object B shell is now a single spherical shell block (replaces 6 wedges)
    block_names_.emplace_back("ObjectBShell");
    add_object_region("ObjectB", "Cube");  // 6 blocks
    first_outer_shell_block_ += 7;
  }
  add_outer_region("Envelope");  // 10 blocks
  first_outer_shell_block_ += 10;
  for (size_t shell = 0; shell < number_of_outer_shells_; shell++) {
    add_outer_region("OuterShell" + std::to_string(shell));  // 10 blocks
  }

  if ((not use_single_block_a_) and (not is_excised_a_)) {
    add_object_interior("ObjectA");  // 1 block
  }
  if ((not use_single_block_b_) and (not is_excised_b_)) {
    add_object_interior("ObjectB");  // 1 block
  }

  ASSERT(block_names_.size() == number_of_blocks_,
         "Number of block names (" << block_names_.size()
                                   << ") doesn't match number of blocks ("
                                   << number_of_blocks_ << ").");

  // Expand initial refinement and number of grid points over all blocks
  const ExpandOverBlocks<std::array<size_t, 3>> expand_over_blocks{
      block_names_, block_groups_};

  try {
    initial_refinement_ = std::visit(expand_over_blocks, initial_refinement);
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialRefinement': " << error.what());
  }
  try {
    initial_number_of_grid_points_ =
        std::visit(expand_over_blocks, initial_number_of_grid_points);
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialGridPoints': " << error.what());
  }

  // For the single spherical shell block replacing the B object shell, the
  // angular dimensions (colatitude = dim 1, longitude = dim 2) use
  // S2Colatitude and S2Longitude topology, which have no boundaries. These
  // dimensions must have refinement level 0; the angular resolution is
  // determined solely by the number of grid points.
  if (not use_single_block_b_) {
    const size_t b_shell_idx = use_single_block_a_ ? 1_st : 12_st;
    initial_refinement_[b_shell_idx][1] = 0;
    initial_refinement_[b_shell_idx][2] = 0;
  }

  // Build time-dependent maps
  std::optional<std::array<double, 3>> radii_A{};
  std::optional<std::array<double, 3>> radii_B{};

  if (not use_single_block_a_) {
    radii_A = std::array{std::get<Object>(object_A_).inner_radius,
                         std::get<Object>(object_A_).outer_radius,
                         sqrt(3.0) * 0.5 * length_inner_cube_};
  }
  if (not use_single_block_b_) {
    radii_B = std::array{std::get<Object>(object_B_).inner_radius,
                         std::get<Object>(object_B_).outer_radius,
                         sqrt(3.0) * 0.5 * length_inner_cube_};
  }

  if (time_dependent_options_.has_value()) {
    // The reason we don't just always use half the inner cube length is to
    // avoid potential roundoff issues if there is no offset
    const std::optional<std::array<double, 3>> cube_A_center =
        length_inner_cube_ == x_coord_a_ - x_coord_b_
            ? std::optional<std::array<double, 3>>{}
            : std::array{translation_ + 0.5 * length_inner_cube_,
                         center_of_mass_offset_[0], center_of_mass_offset_[1]};
    const std::optional<std::array<double, 3>> cube_B_center =
        length_inner_cube_ == x_coord_a_ - x_coord_b_
            ? std::optional<std::array<double, 3>>{}
            : std::array{translation_ - 0.5 * length_inner_cube_,
                         center_of_mass_offset_[0], center_of_mass_offset_[1]};
    time_dependent_options_->build_maps(
        std::array{std::array{x_coord_a_, center_of_mass_offset_[0],
                              center_of_mass_offset_[1]},
                   std::array{x_coord_b_, center_of_mass_offset_[0],
                              center_of_mass_offset_[1]}},
        cube_A_center, cube_B_center,
        std::array{translation_, center_of_mass_offset[0],
                   center_of_mass_offset[1]},
        radii_A, radii_B, not is_excised_a_, not is_excised_b_,
        envelope_radius_, outer_radius_);
  }
}

template <bool UseWorldtube>
Domain<3> BinaryCompactObject<UseWorldtube>::create_domain() const {
  const double inner_sphericity_A = is_excised_a_ ? 1.0 : 0.0;
  // const double inner_sphericity_B = is_excised_b_ ? 1.0 : 0.0;

  using Maps = std::vector<std::unique_ptr<
      CoordinateMapBase<Frame::BlockLogical, Frame::Inertial, 3>>>;

  const std::vector<domain::CoordinateMaps::Distribution>
      object_A_radial_distribution{
          ((not use_single_block_a_) and
           std::get<Object>(object_A_).use_logarithmic_map)
              ? domain::CoordinateMaps::Distribution::Logarithmic
              : domain::CoordinateMaps::Distribution::Linear};

  const std::vector<domain::CoordinateMaps::Distribution>
      object_B_radial_distribution{
          ((not use_single_block_b_) and
           std::get<Object>(object_B_).use_logarithmic_map)
              ? domain::CoordinateMaps::Distribution::Logarithmic
              : domain::CoordinateMaps::Distribution::Linear};

  Maps maps{};

  // ObjectA/B is on the right/left, respectively.
  const Affine3D translation_A{
      Affine{-1.0, 1.0, -1.0 + x_coord_a_ - offset_x_coord_a_,
             1.0 + x_coord_a_ - offset_x_coord_a_},
      Affine{-1.0, 1.0, -1.0 + center_of_mass_offset_[0],
             1.0 + center_of_mass_offset_[0]},
      Affine{-1.0, 1.0, -1.0 + center_of_mass_offset_[1],
             1.0 + center_of_mass_offset_[1]}};
  const Affine3D translation_B{
      Affine{-1.0, 1.0, -1.0 + x_coord_b_ - offset_x_coord_b_,
             1.0 + x_coord_b_ - offset_x_coord_b_},
      Affine{-1.0, 1.0, -1.0 + center_of_mass_offset_[0],
             1.0 + center_of_mass_offset_[0]},
      Affine{-1.0, 1.0, -1.0 + center_of_mass_offset_[1],
             1.0 + center_of_mass_offset_[1]}};

  // Two blocks covering the compact objects and their immediate neighborhood
  if (use_single_block_a_) {
    maps.emplace_back(
        make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(Affine3D{
            Affine(-1.0, 1.0,
                   -0.5 * length_inner_cube_ + x_coord_a_ - offset_x_coord_a_,
                   0.5 * length_inner_cube_ + x_coord_a_ - offset_x_coord_a_),
            Affine(-1.0, 1.0,
                   -0.5 * length_inner_cube_ + center_of_mass_offset_[0],
                   0.5 * length_inner_cube_ + center_of_mass_offset_[0]),
            Affine(-1.0, 1.0,
                   -0.5 * length_inner_cube_ + center_of_mass_offset_[1],
                   0.5 * length_inner_cube_ + center_of_mass_offset_[1])}));
  } else {
    // --- Blocks enclosing each object (12 blocks per object) ---
    //
    // Each object is surrounded by 6 inner wedges that make a sphere, and
    // another 6 outer wedges that transition to a cube.
    const auto& object_a = std::get<Object>(object_A_);
    const auto& offset_a_optional =
        offset_x_coord_a_ == 0
            ? std::nullopt
            : std::make_optional(std::make_pair(
                  length_inner_cube_ * 0.5,
                  std::array<double, 3>{{offset_x_coord_a_, 0.0, 0.0}}));
    Maps maps_center_A =
        domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                Frame::Inertial, 3>(
            sph_wedge_coordinate_maps(
                object_a.inner_radius, object_a.outer_radius,
                inner_sphericity_A, 1.0, use_equiangular_map_,
                offset_a_optional, false, {}, object_A_radial_distribution),
            translation_A);
    Maps maps_cube_A =
        domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                Frame::Inertial, 3>(
            sph_wedge_coordinate_maps(
                object_a.outer_radius, sqrt(3.0) * 0.5 * length_inner_cube_,
                1.0, 0.0, use_equiangular_map_, offset_a_optional),
            translation_A);
    std::move(maps_center_A.begin(), maps_center_A.end(),
              std::back_inserter(maps));
    std::move(maps_cube_A.begin(), maps_cube_A.end(), std::back_inserter(maps));
  }
  if (use_single_block_b_) {
    maps.emplace_back(
        make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(Affine3D{
            Affine(-1.0, 1.0,
                   -0.5 * length_inner_cube_ + x_coord_b_ - offset_x_coord_b_,
                   0.5 * length_inner_cube_ + x_coord_b_ - offset_x_coord_b_),
            Affine(-1.0, 1.0,
                   -0.5 * length_inner_cube_ + center_of_mass_offset_[0],
                   0.5 * length_inner_cube_ + center_of_mass_offset_[0]),
            Affine(-1.0, 1.0,
                   -0.5 * length_inner_cube_ + center_of_mass_offset_[1],
                   0.5 * length_inner_cube_ + center_of_mass_offset_[1])}));
  } else {
    // --- Blocks enclosing object B (1 spherical shell + 6 cube wedges) ---
    //
    // Object B's shell is a single spherical shell block using the
    // SphericalToCartesianPfaffian map (instead of 6 wedge blocks). The cube
    // is still 6 wedge blocks transitioning from sphere to cube.
    const auto& object_b = std::get<Object>(object_B_);
    const auto& offset_b_optional =
        offset_x_coord_b_ == 0
            ? std::nullopt
            : std::make_optional(std::make_pair(
                  length_inner_cube_ * 0.5,
                  std::array<double, 3>{{offset_x_coord_b_, 0.0, 0.0}}));

    // Single spherical shell block for object B shell. Centered directly at
    // the object B position (not offset by cube_scale).
    const Affine3D translation_B_shell{
        Affine{-1.0, 1.0, -1.0 + x_coord_b_, 1.0 + x_coord_b_},
        Affine{-1.0, 1.0, -1.0 + center_of_mass_offset_[0],
               1.0 + center_of_mass_offset_[0]},
        Affine{-1.0, 1.0, -1.0 + center_of_mass_offset_[1],
               1.0 + center_of_mass_offset_[1]}};
    const CoordinateMaps::Interval radial_map_B{-1.0,
                                                1.0,
                                                object_b.inner_radius,
                                                object_b.outer_radius,
                                                object_B_radial_distribution[0],
                                                0.0};
    maps.emplace_back(
        make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
            CoordinateMaps::ProductOf2Maps<CoordinateMaps::Interval,
                                           CoordinateMaps::Identity<2>>{
                radial_map_B, CoordinateMaps::Identity<2>{}},
            CoordinateMaps::SphericalToCartesianPfaffian{},
            translation_B_shell));

    // 6 cube wedge blocks for object B (sphere-to-cube transition)
    Maps maps_cube_B =
        domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                Frame::Inertial, 3>(
            sph_wedge_coordinate_maps(
                object_b.outer_radius, sqrt(3.0) * 0.5 * length_inner_cube_,
                1.0, 0.0, use_equiangular_map_, offset_b_optional),
            translation_B);
    std::move(maps_cube_B.begin(), maps_cube_B.end(), std::back_inserter(maps));
  }

  // --- Frustums enclosing both objects (10 blocks) ---
  //
  // The two abutting cubes are enclosed by a layer of bulged frustums that form
  // a sphere enveloping both objects. While the two objects can be offset from
  // the origin to account for their center of mass, the enveloping frustums are
  // centered at the origin.
  Maps maps_frustums = domain::make_vector_coordinate_map_base<
      Frame::BlockLogical, Frame::Inertial, 3>(frustum_coordinate_maps(
      length_inner_cube_, length_outer_cube_, use_equiangular_map_,
      use_equiangular_map_ and not use_single_block_a_ and
          not use_single_block_b_,
      {{-translation_, -center_of_mass_offset_[0], -center_of_mass_offset_[1]}},
      radial_distribution_envelope_,
      radial_distribution_envelope_ ==
              domain::CoordinateMaps::Distribution::Projective
          ? std::optional<double>(length_inner_cube_ / length_outer_cube_)
          : std::optional<double>(-(length_outer_cube_ + length_inner_cube_) /
                                  (length_outer_cube_ - length_inner_cube_)),
      1.0, opening_angle_));
  std::move(maps_frustums.begin(), maps_frustums.end(),
            std::back_inserter(maps));

  // --- Outer spherical shell (10*num_outer_shells blocks) ---
  Maps maps_outer_shell = domain::make_vector_coordinate_map_base<
      Frame::BlockLogical, Frame::Inertial, 3>(sph_wedge_coordinate_maps(
      envelope_radius_, outer_radius_, 1.0, 1.0, use_equiangular_map_,
      std::nullopt, true, radial_partitioning_outer_shell_,
      radial_distribution_outer_shell_, ShellWedges::All, opening_angle_));
  std::move(maps_outer_shell.begin(), maps_outer_shell.end(),
            std::back_inserter(maps));

  // --- (Optional) object centers (0 to 2 blocks) ---
  //
  // Each object can optionally be filled with a cube-shaped block, in which
  // case the enclosing wedges configured above transition from the cube to a
  // sphere.
  std::unordered_map<std::string, ExcisionSphere<3>> excision_spheres{};
  if (not use_single_block_a_) {
    if (not is_excised_a_) {
      const double scaled_r_inner_A =
          std::get<Object>(object_A_).inner_radius / sqrt(3.0);
      if (use_equiangular_map_) {
        maps.emplace_back(
            make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
                Equiangular3D{
                    Equiangular(-1.0, 1.0,
                                -1.0 * scaled_r_inner_A + offset_x_coord_a_,
                                scaled_r_inner_A + offset_x_coord_a_),
                    Equiangular(-1.0, 1.0, -1.0 * scaled_r_inner_A,
                                scaled_r_inner_A),
                    Equiangular(-1.0, 1.0, -1.0 * scaled_r_inner_A,
                                scaled_r_inner_A)},
                translation_A));
      } else {
        maps.emplace_back(make_coordinate_map_base<Frame::BlockLogical,
                                                   Frame::Inertial>(
            Affine3D{
                Affine(-1.0, 1.0, -1.0 * scaled_r_inner_A + offset_x_coord_a_,
                       scaled_r_inner_A + offset_x_coord_a_),
                Affine(-1.0, 1.0, -1.0 * scaled_r_inner_A, scaled_r_inner_A),
                Affine(-1.0, 1.0, -1.0 * scaled_r_inner_A, scaled_r_inner_A)},
            translation_A));
      }
    }
    // Excision spheres
    // - Block 0 through 5 enclose object A, and 12 through 17 enclose object B.
    // - The 3D wedge map is oriented such that the lower-zeta logical direction
    //   points radially inward.
    else {
      excision_spheres.emplace(
          "ExcisionSphereA",
          ExcisionSphere<3>{std::get<Object>(object_A_).inner_radius,
                            tnsr::I<double, 3, Frame::Grid>{
                                {x_coord_a_, center_of_mass_offset_[0],
                                 center_of_mass_offset_[1]}},
                            {{0, Direction<3>::lower_zeta()},
                             {1, Direction<3>::lower_zeta()},
                             {2, Direction<3>::lower_zeta()},
                             {3, Direction<3>::lower_zeta()},
                             {4, Direction<3>::lower_zeta()},
                             {5, Direction<3>::lower_zeta()}}});
    }
  }
  if (not use_single_block_b_) {
    if (not is_excised_b_) {
      const double scaled_r_inner_B =
          std::get<Object>(object_B_).inner_radius / sqrt(3.0);
      if (use_equiangular_map_) {
        maps.emplace_back(
            make_coordinate_map_base<Frame::BlockLogical, Frame::Inertial>(
                Equiangular3D{
                    Equiangular(-1.0, 1.0,
                                -1.0 * scaled_r_inner_B + offset_x_coord_b_,
                                scaled_r_inner_B + offset_x_coord_b_),
                    Equiangular(-1.0, 1.0, -1.0 * scaled_r_inner_B,
                                scaled_r_inner_B),
                    Equiangular(-1.0, 1.0, -1.0 * scaled_r_inner_B,
                                scaled_r_inner_B)},
                translation_B));
      } else {
        maps.emplace_back(make_coordinate_map_base<Frame::BlockLogical,
                                                   Frame::Inertial>(
            Affine3D{
                Affine(-1.0, 1.0, -1.0 * scaled_r_inner_B + offset_x_coord_b_,
                       scaled_r_inner_B + offset_x_coord_b_),
                Affine(-1.0, 1.0, -1.0 * scaled_r_inner_B, scaled_r_inner_B),
                Affine(-1.0, 1.0, -1.0 * scaled_r_inner_B, scaled_r_inner_B)},
            translation_B));
      }
    }
    // Excision sphere B
    // - Block 0 through 5 enclose object A (wedge lower-zeta points inward).
    // - Object B shell is now a single spherical shell block. The lower-xi
    //   direction of the spherical shell map points radially inward.
    else {
      const size_t first_block_b = use_single_block_a_ ? 1 : 12;
      excision_spheres.emplace(
          "ExcisionSphereB",
          ExcisionSphere<3>{std::get<Object>(object_B_).inner_radius,
                            tnsr::I<double, 3, Frame::Grid>{
                                {x_coord_b_, center_of_mass_offset_[0],
                                 center_of_mass_offset_[1]}},
                            {{first_block_b, Direction<3>::lower_xi()}}});
    }
  }

  // Build block neighbors from the coordinate maps automatically for all
  // conforming (hypercube-topology) interfaces, then manually add the
  // nonconforming interface between the single object B spherical shell block
  // and the 6 object B cube wedge blocks (mirroring
  // NonconformingSphericalShells).
  std::vector<DirectionMap<3, BlockNeighbors<3>>> neighbors_of_all_blocks;
  set_internal_boundaries<3>(make_not_null(&neighbors_of_all_blocks), maps);

  if (not use_single_block_b_) {
    // set_internal_boundaries evaluates all block maps at their logical corners
    // to find conforming neighbors. For the spherical B shell block the map
    // (ProductOf2Maps<Interval,Identity<2>> + SphericalToCartesianPfaffian +
    // translation) can accidentally produce corner positions that coincide with
    // corners of other blocks, causing set_internal_boundaries to add spurious
    // neighbors in directions (e.g. upper_eta) where the spherical_shell
    // topology has no boundary, which would trigger an assert in
    // Element::Element. To prevent this, clear any entries that were set for
    // the B shell block and remove back-pointers from other blocks.
    const size_t b_shell_id = use_single_block_a_ ? 1 : 12;

    // Remove back-pointers from any block that was spuriously connected to the
    // B shell by set_internal_boundaries.
    for (size_t i = 0; i < neighbors_of_all_blocks.size(); ++i) {
      if (i == b_shell_id) {
        continue;
      }
      std::vector<Direction<3>> dirs_to_erase{};
      for (const auto& [dir, nbr] : neighbors_of_all_blocks[i]) {
        if (nbr.ids().contains(b_shell_id)) {
          dirs_to_erase.push_back(dir);
        }
      }
      for (const auto& dir : dirs_to_erase) {
        neighbors_of_all_blocks[i].erase(dir);
      }
    }
    // Clear all spurious neighbors on the B shell itself.
    neighbors_of_all_blocks[b_shell_id].clear();

    // Now set the correct nonconforming neighbors.
    // The B shell's upper_xi face (outer sphere at object_b.outer_radius)
    // is the nonconforming interface with all 6 B cube wedge blocks' lower_zeta
    // faces (inner sphere at the same radius). The orientation map follows
    // NonconformingSphericalShells: the shell's xi (radial) direction maps to
    // the wedge's zeta (radial) direction, with angular dims nonconforming.
    const OrientationMap<3> shell_to_cube{
        {{Direction<3>::upper_zeta(), Direction<3>::self(),
          Direction<3>::self()}}};

    std::unordered_set<size_t> cube_ids{};
    std::unordered_map<size_t, OrientationMap<3>> cube_orientations{};
    for (size_t i = 1; i <= 6; ++i) {
      cube_ids.insert(b_shell_id + i);
      cube_orientations.emplace(b_shell_id + i, shell_to_cube);
    }
    // B shell upper_xi → all 6 cube blocks (nonconforming multi-block)
    neighbors_of_all_blocks[b_shell_id].emplace(
        Direction<3>::upper_xi(),
        BlockNeighbors<3>{cube_ids, cube_orientations, false});
    // Each B cube block lower_zeta → B shell block (nonconforming single)
    const OrientationMap<3> cube_to_shell = shell_to_cube.inverse_map();
    for (size_t i = 1; i <= 6; ++i) {
      neighbors_of_all_blocks[b_shell_id + i].emplace(
          Direction<3>::lower_zeta(),
          BlockNeighbors<3>{
              {b_shell_id}, {{b_shell_id, cube_to_shell}}, false});
    }
  }

  // Build explicit Block objects. The B shell block uses spherical_shell
  // topology so that its angular dimensions (colatitude, longitude) are
  // treated as periodic. All other blocks use the default hypercube topology.
  const size_t b_shell_block_id = use_single_block_a_ ? 1 : 12;
  std::vector<Block<3>> blocks;
  blocks.reserve(number_of_blocks_);
  for (size_t i = 0; i < number_of_blocks_; ++i) {
    const bool is_b_shell =
        (not use_single_block_b_) and (i == b_shell_block_id);
    blocks.emplace_back(std::move(maps[i]), i,
                        std::move(neighbors_of_all_blocks[i]), block_names_[i],
                        is_b_shell ? domain::topologies::spherical_shell
                                   : domain::topologies::hypercube<3>);
  }
  Domain<3> domain{std::move(blocks), std::move(excision_spheres),
                   block_groups_};

  // Inject the hard-coded time-dependence
  if (time_dependent_options_.has_value()) {
    // Default initialize everything to nullptr so that we only need to set the
    // appropriate block maps for the specific frames
    std::vector<std::unique_ptr<
        domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>>>
        grid_to_inertial_block_maps{number_of_blocks_};
    std::vector<std::unique_ptr<
        domain::CoordinateMapBase<Frame::Grid, Frame::Distorted, 3>>>
        grid_to_distorted_block_maps{number_of_blocks_};
    std::vector<std::unique_ptr<
        domain::CoordinateMapBase<Frame::Distorted, Frame::Inertial, 3>>>
        distorted_to_inertial_block_maps{number_of_blocks_};

    // Some maps (e.g. expansion, rotation) are applied to all blocks,
    // while other maps (e.g. size) are only applied to some blocks. Also, some
    // maps are applied from the Grid to the Inertial frame, while others are
    // applied to/from the Distorted frame. So there are several different
    // distinct combinations of time-dependent maps that will be applied. This
    // should largely be taken care of by the TimeDependentOptions.

    // Single blocks are added after the outer shell if an object's interior is
    // not excised, so to find the final outer shell block, subtract the single
    // blocks that may have been added.
    size_t final_block_outer_shell = number_of_blocks_ - 1;
    if ((not use_single_block_a_) and (not is_excised_a_)) {
      --final_block_outer_shell;
    }
    if ((not use_single_block_b_) and (not is_excised_b_)) {
      --final_block_outer_shell;
    }

    // All blocks except possibly the first 6 or 12 blocks of each object get
    // the same map from the Grid to the Inertial frame, so initialize the final
    // block with the "base" map (here a composition of an expansion and a
    // rotation). When covering the inner regions with cubes, all blocks will
    // use the same time-dependent map instead.
    grid_to_inertial_block_maps[final_block_outer_shell] =
        time_dependent_options_
            ->grid_to_inertial_map<domain::ObjectLabel::None>(std::nullopt,
                                                              false);
    if (final_block_outer_shell < number_of_blocks_ - 1) {
      grid_to_inertial_block_maps[number_of_blocks_ - 1] =
          time_dependent_options_
              ->grid_to_inertial_map<domain::ObjectLabel::None>(std::nullopt,
                                                                true);
    }
    size_t final_block_envelope = first_outer_shell_block_ - 1;

    grid_to_inertial_block_maps[final_block_envelope] =
        time_dependent_options_
            ->grid_to_inertial_map<domain::ObjectLabel::None>(std::nullopt,
                                                              true);
    // Inside the excision sphere we add a special grid to inertial map that can
    // be evaluated at the center. This allows the center of the
    // excisions/horizons to be mapped properly to the inertial frame.
    if (is_excised_a_ and
        grid_to_inertial_block_maps[number_of_blocks_ - 1] != nullptr) {
      domain.inject_time_dependent_map_for_excision_sphere(
          "ExcisionSphereA",
          time_dependent_options_->grid_to_inertial_map<domain::ObjectLabel::A>(
              {0}, true, true));
    }
    if (is_excised_b_ and
        grid_to_inertial_block_maps[number_of_blocks_ - 1] != nullptr) {
      domain.inject_time_dependent_map_for_excision_sphere(
          "ExcisionSphereB",
          time_dependent_options_->grid_to_inertial_map<domain::ObjectLabel::B>(
              {0}, true, true));
    }

    const size_t first_block_object_B = use_single_block_a_ ? 1 : 12;

    // We loop over all blocks. If we are using a single block for either A or
    // B, then we only need a grid to inertial map; there is no distorted frame.
    // If we don't have a single block around A, then for 12 blocks around
    // object A we determine if there is a central cube. If there is, then
    // there is no distorted frame (this is signaled by
    // block_for_distorted_frame being a nullopt). If the region is excised,
    // then it is up to the time dependent options to know if a distorted frame
    // exists in that block. Therefore we just pass the relative block number to
    // the time_dependent_options_ functions. The remaining blocks only get the
    // grid to inertial map.
    //
    // Object B shell is now a single spherical shell block (block index 0
    // relative to B), followed by 6 cube wedge blocks (relative indices 1-6).
    // The time-dependent options expect relative block indices 0-5 for shell
    // blocks (with shape map) and 6-11 for cube blocks (no shape map). We map:
    //   - Spherical shell (relative 0) -> index 0 (gets shape map since < 6)
    //   - Cube blocks (relative 1-6)   -> indices 6-11 (no shape map since >=
    //   6)
    for (size_t block = 0; block < number_of_blocks_ - 1; ++block) {
      if ((not use_single_block_a_) and block < first_block_object_B) {
        const size_t block_for_distorted_frame = block;
        grid_to_inertial_block_maps[block] =
            time_dependent_options_
                ->grid_to_inertial_map<domain::ObjectLabel::A>(
                    block_for_distorted_frame, true);
        grid_to_distorted_block_maps[block] =
            time_dependent_options_
                ->grid_to_distorted_map<domain::ObjectLabel::A>(
                    block_for_distorted_frame);
        distorted_to_inertial_block_maps[block] =
            time_dependent_options_
                ->distorted_to_inertial_map<domain::ObjectLabel::A>(
                    block_for_distorted_frame, true);
      } else if ((not use_single_block_b_) and block >= first_block_object_B and
                 block < first_block_object_B + 7) {
        // Spherical shell block (relative 0) gets index 0.
        // Cube blocks (relative 1-6) get indices 6-11 to avoid shape map.
        const size_t block_relative = block - first_block_object_B;
        const size_t block_for_distorted_frame =
            block_relative == 0 ? 0 : block_relative + 5;
        grid_to_inertial_block_maps[block] =
            time_dependent_options_
                ->grid_to_inertial_map<domain::ObjectLabel::B>(
                    block_for_distorted_frame, true);
        grid_to_distorted_block_maps[block] =
            time_dependent_options_
                ->grid_to_distorted_map<domain::ObjectLabel::B>(
                    block_for_distorted_frame);
        distorted_to_inertial_block_maps[block] =
            time_dependent_options_
                ->distorted_to_inertial_map<domain::ObjectLabel::B>(
                    block_for_distorted_frame, true);
        // check if block is less than outer shell block for rigid
        // expansion/translation, but no distorted map
      } else if (block < first_outer_shell_block_ and
                 grid_to_inertial_block_maps[number_of_blocks_ - 1] !=
                     nullptr) {
        grid_to_inertial_block_maps[block] =
            grid_to_inertial_block_maps[final_block_envelope]->get_clone();
      } else if (block > final_block_outer_shell and
                 grid_to_inertial_block_maps[number_of_blocks_ - 1] !=
                     nullptr) {
        // the inner cube blocks are after outershell and we want to copy the
        // corresponding object blocks.
        if ((not use_single_block_a_) and (not is_excised_a_)) {
          grid_to_inertial_block_maps[block] =
              grid_to_inertial_block_maps[0]->get_clone();
        }
        if ((not use_single_block_b_) and (not is_excised_b_)) {
          grid_to_inertial_block_maps[block] =
              grid_to_inertial_block_maps[first_block_object_B]->get_clone();
        }
      } else if (grid_to_inertial_block_maps[number_of_blocks_ - 1] !=
                 nullptr) {
        // No distorted frame
        grid_to_inertial_block_maps[block] =
            grid_to_inertial_block_maps[final_block_outer_shell]->get_clone();
      }
    }
    // Finally, inject the time dependent maps into the corresponding blocks
    for (size_t block = 0; block < number_of_blocks_; ++block) {
      if (grid_to_inertial_block_maps[block] == nullptr) {
        continue;
      }
      domain.inject_time_dependent_map_for_block(
          block, std::move(grid_to_inertial_block_maps[block]),
          std::move(grid_to_distorted_block_maps[block]),
          std::move(distorted_to_inertial_block_maps[block]));
    }
  }
  return domain;
}

template <bool UseWorldtube>
std::vector<DirectionMap<
    3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
BinaryCompactObject<UseWorldtube>::external_boundary_conditions() const {
  if (outer_boundary_condition_ == nullptr) {
    return {};
  }
  std::vector<DirectionMap<
      3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
      boundary_conditions{number_of_blocks_};
  // Excision surfaces
  // Blocks 0-5 are object A's 6 shell wedges; lower-zeta points radially
  // inward toward the excision.
  if (is_excised_a_) {
    for (size_t i = 0; i < 6; ++i) {
      boundary_conditions[i][Direction<3>::lower_zeta()] =
          (*(std::get<Object>(object_A_).inner_boundary_condition))
              ->get_clone();
    }
  }
  // Object B's shell is now a single spherical shell block; lower-xi points
  // radially inward toward the excision.
  const size_t first_block_object_B = use_single_block_a_ ? 1 : 12;
  if (is_excised_b_) {
    boundary_conditions[first_block_object_B][Direction<3>::lower_xi()] =
        (*(std::get<Object>(object_B_).inner_boundary_condition))->get_clone();
  }
  // Outer boundary
  // Block layout (without single-block replacements):
  //   0-11 : Object A shell+cube (12 blocks)
  //   12   : Object B spherical shell (1 block)
  //   13-18: Object B cube (6 blocks)
  //   19-28: Envelope (10 blocks)
  //   29+  : Outer shells
  const size_t offset_outer_blocks =
      (use_single_block_a_ and use_single_block_b_
           ? 12
           : (use_single_block_a_ ? 18 : (use_single_block_b_ ? 23 : 29))) +
      10 * (number_of_outer_shells_ - 1);
  for (size_t i = 0; i < 10; ++i) {
    boundary_conditions[i + offset_outer_blocks][Direction<3>::upper_zeta()] =
        outer_boundary_condition_->get_clone();
  }
  return boundary_conditions;
}

template <bool UseWorldtube>
std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
BinaryCompactObject<UseWorldtube>::functions_of_time(
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const {
  return time_dependent_options_.has_value()
             ? time_dependent_options_->create_functions_of_time<UseWorldtube>(
                   initial_expiration_times)
             : std::unordered_map<
                   std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>{};
}

template class BinaryCompactObject<true>;
template class BinaryCompactObject<false>;
}  // namespace domain::creators
