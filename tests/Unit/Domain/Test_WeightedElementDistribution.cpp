// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Domain/Block.hpp"
#include "Domain/Creators/AlignedLattice.hpp"
#include "Domain/Creators/BinaryCompactObject.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Domain.hpp"
#include "Domain/OptionTags.hpp"
#include "Domain/Protocols/Metavariables.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/WeightedElementDistribution.hpp"
#include "Domain/ZCurve.hpp"
#include "Framework/TestCreation.hpp"
#include "Helpers/Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {
// Make a `domain::creators::AlignedLattice` from an option string
template <size_t Dim>
auto make_domain_creator(const std::string& opt_string) {
  return TestHelpers::test_option_tag<
      domain::OptionTags::DomainCreator<Dim>,
      TestHelpers::domain::BoundaryConditions::
          MetavariablesWithoutBoundaryConditions<
              Dim, domain::creators::AlignedLattice<Dim>>>(opt_string);
}

// Metavariables for a `domain::creators::BinaryCompactObject`
template <size_t Dim, bool EnableTimeDependentMaps, bool WithBoundaryConditions>
struct Metavariables {
  struct domain : tt::ConformsTo<::domain::protocols::Metavariables> {
    static constexpr bool enable_time_dependent_maps = EnableTimeDependentMaps;
  };
  using system = tmpl::conditional_t<WithBoundaryConditions,
                                     TestHelpers::domain::BoundaryConditions::
                                         SystemWithBoundaryConditions<Dim>,
                                     TestHelpers::domain::BoundaryConditions::
                                         SystemWithoutBoundaryConditions<Dim>>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<tmpl::pair<
        DomainCreator<3>, tmpl::list<::domain::creators::BinaryCompactObject>>>;
  };
};

// Stringify a bool
std::string stringize(const bool t) { return t ? "true" : "false"; }

// Create an option string for a `domain::creators::BinaryCompactObject`
std::string create_option_string(const bool excise_A, const bool excise_B,
                                 const bool add_time_dependence,
                                 const bool use_logarithmic_map_AB,
                                 const bool use_equiangular_map,
                                 const size_t additional_refinement_outer,
                                 const size_t additional_refinement_A,
                                 const size_t additional_refinement_B,
                                 const bool add_boundary_condition) {
  const std::string time_dependence{
      add_time_dependence ? "  TimeDependentMaps:\n"
                            "    InitialTime: 1.0\n"
                            "    ExpansionMap: \n"
                            "      OuterBoundary: 25.0\n"
                            "      InitialExpansion: 1.0\n"
                            "      InitialExpansionVelocity: -0.1\n"
                            "      AsymptoticVelocityOuterBoundary: -0.1\n"
                            "      DecayTimescaleOuterBoundaryVelocity: 5.0\n"
                            "    RotationMap:\n"
                            "      InitialAngularVelocity: [0.0, 0.0, -0.2]\n"
                            "    SizeMap:\n"
                            "      InitialValues: [0.0, 0.0]\n"
                            "      InitialVelocities: [-0.1, -0.2]\n"
                            "      InitialAccelerations: [0.01, 0.02]"
                          : ""};
  const std::string interior_A{
      add_boundary_condition
          ? std::string{"    Interior:\n" +
                        std::string{excise_A
                                        ? "      ExciseWithBoundaryCondition:\n"
                                          "        TestBoundaryCondition:\n"
                                          "          Direction: lower-zeta\n"
                                          "          BlockId: 50\n"
                                        : "      Auto\n"}}
          : "    ExciseInterior: " + stringize(excise_A) + "\n"};
  const std::string interior_B{
      add_boundary_condition
          ? std::string{"    Interior:\n" +
                        std::string{excise_B
                                        ? "      ExciseWithBoundaryCondition:\n"
                                          "        TestBoundaryCondition:\n"
                                          "          Direction: lower-zeta\n"
                                          "          BlockId: 50\n"
                                        : "      Auto\n"}}
          : "    ExciseInterior: " + stringize(excise_B) + "\n"};
  const std::string outer_boundary_condition{
      add_boundary_condition ? std::string{"    BoundaryCondition:\n"
                                           "      TestBoundaryCondition:\n"
                                           "        Direction: upper-zeta\n"
                                           "        BlockId: 50\n"}
                             : ""};
  return "BinaryCompactObject:\n"
         "  ObjectA:\n"
         "    InnerRadius: 1.0\n"
         "    OuterRadius: 2.0\n"
         "    XCoord: 3.0\n" +
         interior_A +
         "    UseLogarithmicMap: " + stringize(use_logarithmic_map_AB) +
         "\n"
         "  ObjectB:\n"
         "    InnerRadius: 0.2\n"
         "    OuterRadius: 1.0\n"
         "    XCoord: -2.0\n" +
         interior_B +
         "    UseLogarithmicMap: " + stringize(use_logarithmic_map_AB) +
         "\n"
         "  EnvelopingCube:\n"
         "    Radius: 22.0\n"
         "    UseProjectiveMap: true\n"
         "    Sphericity: 1.0\n"
         "  OuterShell:\n"
         "    InnerRadius: Auto\n"
         "    OuterRadius: 25.0\n"
         "    RadialDistribution: Linear\n" +
         outer_boundary_condition + "  InitialRefinement:\n" +
         (excise_A ? "" : "    ObjectAInterior: [1, 1, 1]\n") +
         (excise_B ? "" : "    ObjectBInterior: [1, 1, 1]\n") +
         "    ObjectAShell: [1, 1, " +
         std::to_string(1 + additional_refinement_A) +
         "]\n"
         "    ObjectBShell: [1, 1, " +
         std::to_string(1 + additional_refinement_B) +
         "]\n"
         "    ObjectACube: [1, 1, 1]\n"
         "    ObjectBCube: [1, 1, 1]\n"
         "    EnvelopingCube: [1, 1, 1]\n"
         "    OuterShell: [1, 1, " +
         std::to_string(1 + additional_refinement_outer) +
         "]\n"
         "  InitialGridPoints: 3\n"
         "  UseEquiangularMap: " +
         stringize(use_equiangular_map) + "\n" + time_dependence;
}

// Test the computation of the weighting done by
// `domain::WeightedBlockZCurveProcDistribution::get_cost_by_element_by_block`
void test_cost_function() {
  const auto domain_creator1 = make_domain_creator<3>(
      "AlignedLattice:\n"
      "  BlockBounds: [[0.0, 1.0, 2.0], [0.0, 1.0], [0.0, 1.0]]\n" +
      std::string{"  IsPeriodicIn: [false, false, false]\n"} +
      "  InitialGridPoints: [4, 4, 4]\n"
      "  InitialLevels: [2, 1, 0]\n"
      "  RefinedLevels: []\n"
      "  RefinedGridPoints: []\n"
      "  BlocksToExclude: []\n");
  const auto* aligned_blocks_creator1 =
      dynamic_cast<const domain::creators::AlignedLattice<3>*>(
          domain_creator1.get());
  const auto domain1 = aligned_blocks_creator1->create_domain();
  const auto& blocks1 = domain1.blocks();

  // Block size and grid points are the same as blocks in `domain1`, but
  // refinement levels are different
  const auto domain_creator2 = make_domain_creator<3>(
      "AlignedLattice:\n"
      "  BlockBounds: [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]\n" +
      std::string{"  IsPeriodicIn: [false, false, false]\n"} +
      "  InitialGridPoints: [4, 4, 4]\n"
      "  InitialLevels: [2, 3, 2]\n"
      "  RefinedLevels: []\n"
      "  RefinedGridPoints: []\n"
      "  BlocksToExclude: []\n");
  const auto* aligned_blocks_creator2 =
      dynamic_cast<const domain::creators::AlignedLattice<3>*>(
          domain_creator2.get());
  const auto domain2 = aligned_blocks_creator2->create_domain();
  const auto& blocks2 = domain2.blocks();

  // Block size and refinement levels are the same as blocks in `domain1`, but
  // grid points are different
  const auto domain_creator3 = make_domain_creator<3>(
      "AlignedLattice:\n"
      "  BlockBounds: [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]\n" +
      std::string{"  IsPeriodicIn: [false, false, false]\n"} +
      "  InitialGridPoints: [4, 3, 2]\n"
      "  InitialLevels: [2, 1, 0]\n"
      "  RefinedLevels: []\n"
      "  RefinedGridPoints: []\n"
      "  BlocksToExclude: []\n");
  const auto* aligned_blocks_creator3 =
      dynamic_cast<const domain::creators::AlignedLattice<3>*>(
          domain_creator3.get());
  const auto domain3 = aligned_blocks_creator3->create_domain();
  const auto& blocks3 = domain2.blocks();

  const auto costs1 = domain::WeightedBlockZCurveProcDistribution<3>::
      get_cost_by_element_by_block(
          blocks1, aligned_blocks_creator1->initial_refinement_levels(),
          aligned_blocks_creator1->initial_extents(),
          Spectral::Quadrature::GaussLobatto);

  const auto costs2 = domain::WeightedBlockZCurveProcDistribution<3>::
      get_cost_by_element_by_block(
          blocks2, aligned_blocks_creator2->initial_refinement_levels(),
          aligned_blocks_creator2->initial_extents(),
          Spectral::Quadrature::GaussLobatto);

  const auto costs3 = domain::WeightedBlockZCurveProcDistribution<3>::
      get_cost_by_element_by_block(
          blocks3, aligned_blocks_creator3->initial_refinement_levels(),
          aligned_blocks_creator3->initial_extents(),
          Spectral::Quadrature::GaussLobatto);

  Approx custom_approx_e16 = Approx::custom().epsilon(1.0e-15).scale(1.0);

  const double elemental_cost1 = costs1[0][0];
  for (size_t i = 0; i < costs1[0].size(); i++) {
    // check that all elements in a block have the same cost
    CHECK(elemental_cost1 == custom_approx_e16(costs1[0][i]));
    CHECK(elemental_cost1 == custom_approx_e16(costs1[1][i]));
  }
  // check that the elements in both blocks have the same cost
  CHECK_ITERABLE_APPROX(costs1[0], costs1[1]);

  const double elemental_cost2 = costs2[0][0];
  for (size_t i = 1; i < costs2[0].size(); i++) {
    // check that all elements in the block have the same cost
    CHECK(elemental_cost2 == costs2[0][i]);
  }

  Approx custom_approx_e15 = Approx::custom().epsilon(1.0e-15).scale(1.0);
  // The highest refinement for the first test domain is 2 while the highest
  // refinement for the second test domain is 3, grid points held constant.
  // Since the minimum grid spacing of the second domain is half the minimum
  // grid spacing of the first and since
  // elemental cost = (# of grid points) / sqrt(min grid spacing), the
  // elemental cost of the second domain should be a factor of sqrt(2) the cost.
  CHECK(elemental_cost2 == custom_approx_e15(sqrt(2.0) * elemental_cost1));

  const double elemental_cost3 = costs3[0][0];
  for (size_t i = 1; i < costs3[0].size(); i++) {
    // check that all elements in the block have the same cost
    CHECK(elemental_cost3 == custom_approx_e16(costs3[0][i]));
  }

  // The minimum grid spacing for the first and third domain are equal, but the
  // number of grid points in an element in the first is 64 while the number of
  // grid points in an element in the third is 24. Since
  // elemental cost = (# of grid points) / sqrt(min grid spacing), the
  // elemental cost of the second domain should be a factor of 24/64 = 3/8 the
  // cost.
  CHECK(elemental_cost3 == elemental_cost1 * 3.0 / 8.0);
}

// Test the processor distribution logic of the
// `domain::WeightedBlockZCurveProcDistribution` constructor
template <size_t Dim>
void test_element_distribution(
    const DomainCreator<Dim>& domain_creator,
    const size_t number_of_procs_with_elements,
    const std::unordered_set<size_t>& global_procs_to_ignore = {}) {
  const auto domain = domain_creator.create_domain();
  const auto& blocks = domain.blocks();
  const auto initial_refinement_levels =
      domain_creator.initial_refinement_levels();
  const auto initial_extents = domain_creator.initial_extents();

  const size_t num_blocks = blocks.size();
  size_t num_elements = 0;
  std::vector<size_t> num_elements_by_block(num_blocks);
  std::fill(num_elements_by_block.begin(), num_elements_by_block.end(), 0);
  for (size_t i = 0; i < num_blocks; i++) {
    size_t num_elements_this_block =
        two_to_the(gsl::at(initial_refinement_levels[i], 0));
    for (size_t j = 1; j < Dim; j++) {
      num_elements_this_block *=
          two_to_the(gsl::at(initial_refinement_levels[i], j));
    }
    num_elements_by_block[i] = num_elements_this_block;
    num_elements += num_elements_this_block;
  }

  const domain::WeightedBlockZCurveProcDistribution<Dim> element_distribution(
      number_of_procs_with_elements, blocks, initial_refinement_levels,
      initial_extents, Spectral::Quadrature::GaussLobatto,
      global_procs_to_ignore);
  const auto proc_map = element_distribution.block_element_distribution();

  const size_t total_procs =
      number_of_procs_with_elements + global_procs_to_ignore.size();
  std::vector<size_t> num_elements_by_proc(total_procs);
  std::fill(num_elements_by_proc.begin(), num_elements_by_proc.end(), 0);
  std::vector<size_t> actual_num_elements_by_block_in_dist(num_blocks);
  std::fill(actual_num_elements_by_block_in_dist.begin(),
            actual_num_elements_by_block_in_dist.end(), 0);
  size_t actual_num_elements_in_dist = 0;
  for (size_t block_number = 0; block_number < proc_map.size();
       block_number++) {
    for (const auto& proc_allowance : proc_map[block_number]) {
      const size_t proc_number = proc_allowance.first;
      const size_t element_allowance = proc_allowance.second;
      num_elements_by_proc[proc_number] += element_allowance;
      actual_num_elements_by_block_in_dist[block_number] += element_allowance;
    }
    // check that the number of elements in this block accounted for in the
    // distribution matches the expected number of total elements for this block
    CHECK(actual_num_elements_by_block_in_dist[block_number] ==
          num_elements_by_block[block_number]);
    actual_num_elements_in_dist +=
        actual_num_elements_by_block_in_dist[block_number];
  }
  // check that the number of elements accounted for in the distribution matches
  // the expected number of total elements
  CHECK(actual_num_elements_in_dist == num_elements);

  const auto costs = domain::WeightedBlockZCurveProcDistribution<
      Dim>::get_cost_by_element_by_block(blocks, initial_refinement_levels,
                                         initial_extents,
                                         Spectral::Quadrature::GaussLobatto);

  double total_cost = 0.0;
  for (const auto& block : costs) {
    for (const double element_cost : block) {
      total_cost += element_cost;
    }
  }

  // one flattened vector instead of vectors by Block
  std::vector<double> costs_flattened(num_elements);
  size_t cost_index = 0;
  for (const auto& block : costs) {
    for (const double element_cost : block) {
      costs_flattened[cost_index] = element_cost;
      cost_index++;
    }
  }

  cost_index = 0;
  double cost_remaining = total_cost;
  size_t procs_skipped = 0;
  // check that we distributed the right number of elements to each proc based
  // on the sum of their costs in Z-curve index order
  for (size_t i = 0; i < total_procs; i++) {
    if (global_procs_to_ignore.count(i)) {
      procs_skipped++;
      continue;
    }

    if (cost_index < num_elements) {
      // if we haven't accounted for all elements yet, we should still have cost
      // left to account for
      CHECK(cost_remaining <= total_cost);
    } else {
      // if we've already accounted for all the elements, we shouldn't have any
      // cost left to account for, and it's the case that more procs were
      // requested than could be used, e.g. in the case of less elements than
      // procs
      Approx custom_approx = Approx::custom().epsilon(1.0e-11).scale(1.0);
      CHECK(cost_remaining == custom_approx(0.0));
      break;
    }

    // the average cost per proc that we're aiming for
    const double target_proc_cost =
        cost_remaining / (number_of_procs_with_elements - i + procs_skipped);

    // total cost on the processor before adding the cost of the final element
    // assigned to this proc
    double proc_cost_without_final_element = 0.0;
    const size_t num_elements_this_proc = num_elements_by_proc[i];
    // add up costs of all elements but the final one to add
    for (size_t j = 0; j < num_elements_this_proc - 1; j++) {
      const double this_cost = costs_flattened[cost_index + j];
      proc_cost_without_final_element += this_cost;
    }

    // the cost of all of the elements assigned to this proc
    const double proc_cost_with_final_element =
        proc_cost_without_final_element +
        costs_flattened[cost_index + num_elements_this_proc - 1];

    const double diff_without_final_element =
        abs(proc_cost_without_final_element - target_proc_cost);

    const double diff_with_final_element =
        abs(proc_cost_with_final_element - target_proc_cost);

    // if the elements assigned to this proc have a cost that is over the target
    // cost per proc, make sure that either it's because only one element is
    // being assigned to the proc or this cost is closer to the target than if
    // we omitted the final element, i.e. check that it's better to keep the
    // final element than to not
    if (proc_cost_with_final_element > target_proc_cost) {
      const bool result = num_elements_this_proc == 1 or
                          diff_with_final_element <= diff_without_final_element;
      CHECK(result);
    }

    if (cost_index + num_elements_this_proc < num_elements) {
      // total cost on the processor if we were to add the cost of the next
      // element (one additional than the number chosen)
      const double proc_cost_with_extra_element =
          proc_cost_with_final_element +
          costs_flattened[cost_index + num_elements_this_proc];
      const double diff_with_extra_element =
          abs(proc_cost_with_extra_element - target_proc_cost);

      // check that it's better to not add one more element than the number
      // chosen
      CHECK(diff_with_extra_element >= diff_with_final_element);
    }

    cost_index += num_elements_this_proc;
    cost_remaining -= proc_cost_with_final_element;
  }

  // check that any remainder of processors we didn't need do indeed have 0
  // elements assigned to them
  for (size_t j = cost_index + 1; j < total_procs; j++) {
    CHECK(num_elements_by_proc[j] == 0);
  }
}

// Test the retrieval of the assigned processor that is done by
// `domain::WeightedBlockZCurveProcDistribution::get_proc_for_element`
template <size_t Dim>
void test_proc_retrieval(
    const DomainCreator<Dim>& domain_creator,
    const size_t number_of_procs_with_elements,
    const std::unordered_set<size_t>& global_procs_to_ignore = {}) {
  const auto domain = domain_creator.create_domain();
  const auto& blocks = domain.blocks();
  const auto initial_refinement_levels =
      domain_creator.initial_refinement_levels();
  const auto initial_extents = domain_creator.initial_extents();

  const size_t num_blocks = blocks.size();
  std::vector<std::vector<ElementId<Dim>>> element_ids_in_z_curve_order(
      num_blocks);
  for (size_t i = 0; i < num_blocks; i++) {
    element_ids_in_z_curve_order[i] =
        domain::initial_element_ids_in_z_curve_order(
            i, gsl::at(initial_refinement_levels, i), 0);
  }

  const domain::WeightedBlockZCurveProcDistribution<Dim> element_distribution(
      number_of_procs_with_elements, blocks, initial_refinement_levels,
      initial_extents, Spectral::Quadrature::GaussLobatto,
      global_procs_to_ignore);
  const auto proc_map = element_distribution.block_element_distribution();

  const size_t total_number_of_procs =
      number_of_procs_with_elements + global_procs_to_ignore.size();

  // whether or not we've assigned elements to a proc
  std::vector<bool> proc_hit(total_number_of_procs);
  std::fill(proc_hit.begin(), proc_hit.end(), false);

  size_t highest_proc_assigned = 0;
  for (size_t i = 0; i < num_blocks; i++) {
    size_t element_index = 0;
    const std::vector<std::pair<size_t, size_t>>& proc_map_this_block =
        proc_map[i];
    const size_t num_procs_this_block = proc_map_this_block.size();

    for (size_t j = 0; j < num_procs_this_block; j++) {
      const size_t expected_proc = proc_map_this_block[j].first;
      const size_t proc_allowance = proc_map_this_block[j].second;

      for (size_t k = 0; k < proc_allowance; k++) {
        const size_t actual_proc = element_distribution.get_proc_for_element(
            element_ids_in_z_curve_order[i][element_index]);
        // check that the correct processor is returned for the `ElementId`
        CHECK(actual_proc == expected_proc);
        proc_hit[actual_proc] = true;
        if (highest_proc_assigned < actual_proc) {
          highest_proc_assigned = actual_proc;
        }
      }
      element_index += proc_allowance;
    }
  }

  // check that ignored procs were indeed skipped and that all other procs
  // up to the highest one assigned were hit
  for (size_t i = 0; i < highest_proc_assigned + 1; i++) {
    if (global_procs_to_ignore.count(i) == 0) {
      CHECK(proc_hit[i]);
    } else {
      CHECK(not proc_hit[i]);
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.WeightedElementDistribution", "[Domain][Unit]") {
  // Test computation of elemental weights
  test_cost_function();

  // Test inputs

  // 1D, single block
  const auto lattice_1d = make_domain_creator<1>(
      "AlignedLattice:\n"
      "  BlockBounds: [[0.0, 1.0]]\n" +
      std::string{"  IsPeriodicIn: [false]\n"} +
      "  InitialGridPoints: [6]\n"
      "  InitialLevels: [4]\n"
      "  RefinedLevels: []\n"
      "  RefinedGridPoints: []\n"
      "  BlocksToExclude: []\n");

  // 2D
  const auto lattice_2d = make_domain_creator<2>(
      "AlignedLattice:\n"
      "  BlockBounds: [[0.0, 0.3], [0.0, 0.8, 2.5, 4.9]]\n" +
      std::string{"  IsPeriodicIn: [false, false]\n"} +
      "  InitialGridPoints: [4, 5]\n"
      "  InitialLevels: [2, 3]\n"
      "  RefinedLevels: []\n"
      "  RefinedGridPoints: []\n"
      "  BlocksToExclude: []\n");

  // 3D
  const auto binary_compact_object_creator =
      TestHelpers::test_option_tag<domain::OptionTags::DomainCreator<3>,
                                   Metavariables<3, true, false>>(
          create_option_string(true, true, true, false, false, 0, 0, 0, false));

  // Test element distribution and proc retrieval for 1D, 2D, and 3D. For each
  // dimension, four cases are tested: single proc requested, multiple procs
  // requested, procs to ignore requested, and more procs requested than
  // elements to distribute.

  // 1D
  test_element_distribution(*lattice_1d, 1);
  test_element_distribution(*lattice_1d, 5);
  test_element_distribution(*lattice_1d, 10, std::unordered_set<size_t>{4, 6});
  test_element_distribution(*lattice_1d, 33, std::unordered_set<size_t>{7});

  // 2D
  test_element_distribution(*lattice_2d, 1);
  test_element_distribution(*lattice_2d, 5);
  test_element_distribution(*lattice_2d, 20, std::unordered_set<size_t>{4, 20});
  test_element_distribution(*lattice_2d, 54, std::unordered_set<size_t>{0, 1});

  // 3D
  test_element_distribution(*binary_compact_object_creator, 1);
  test_element_distribution(*binary_compact_object_creator, 12);
  test_element_distribution(*binary_compact_object_creator, 73,
                            std::unordered_set<size_t>{5, 8, 9, 75});
  test_element_distribution(*binary_compact_object_creator, 500,
                            std::unordered_set<size_t>{100});

  // 1D
  test_proc_retrieval(*lattice_1d, 1);
  test_proc_retrieval(*lattice_1d, 5);
  test_proc_retrieval(*lattice_1d, 10, std::unordered_set<size_t>{4, 6});
  test_proc_retrieval(*lattice_1d, 33, std::unordered_set<size_t>{7});

  // 2D
  test_proc_retrieval(*lattice_2d, 1);
  test_proc_retrieval(*lattice_2d, 5);
  test_proc_retrieval(*lattice_2d, 20, std::unordered_set<size_t>{4, 20});
  test_proc_retrieval(*lattice_2d, 54, std::unordered_set<size_t>{0, 1});

  // 3D
  test_proc_retrieval(*binary_compact_object_creator, 1);
  test_proc_retrieval(*binary_compact_object_creator, 12);
  test_proc_retrieval(*binary_compact_object_creator, 73,
                      std::unordered_set<size_t>{5, 8, 9, 75});
  test_proc_retrieval(*binary_compact_object_creator, 500,
                      std::unordered_set<size_t>{100});
}
