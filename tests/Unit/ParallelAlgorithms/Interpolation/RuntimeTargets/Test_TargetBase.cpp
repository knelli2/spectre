// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>
#include <pup.h>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/AlignedLattice.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Framework/TestHelpers.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Target.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace intrp::Targets {
namespace {
template <size_t Dim>
class TestTarget : public Target<Dim> {
 public:
  TestTarget() = default;

  // Just make it easy and make all the components the same value
  explicit TestTarget(const size_t number_of_grid_points,
                      const double fill_value)
      : points_(number_of_grid_points, fill_value) {}

  explicit TestTarget(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_base_template(Target<Dim>, TestTarget<Dim>);  // NOLINT

  const std::string& name() const override { return name_; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    Target<Dim>::pup(p);
    p | name_;
    p | points_;
  }

 private:
  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
      const override {
    return points_;
  }

  std::string name_{"TestTarget"};
  tnsr::I<DataVector, Dim, Frame::NoFrame> points_{};
};

template <size_t Dim>
PUP::able::PUP_ID TestTarget<Dim>::my_PUP_ID = 0;  // NOLINT

template <size_t Dim>
class TestTargetOperator : public Target<Dim> {
 public:
  TestTargetOperator() = default;

  // Just make it easy and make all the components the same value
  explicit TestTargetOperator(const size_t number_of_grid_points,
                              const double fill_value)
      : points_(number_of_grid_points, fill_value) {}

  explicit TestTargetOperator(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_base_template(Target<Dim>,
                                     TestTargetOperator<Dim>);  // NOLINT

  using argument_tags = tmpl::list<>;

  // We purposefully error so we know the corrector operator() was called
  template <typename Metavariables>
  std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>
  operator()(const Parallel::GlobalCache<Metavariables>& /*cache*/,
             const double /*time*/, const std::string& /*frame*/) const {
    ERROR("Called the correct operator()");
  }

  const std::string& name() const override { return name_; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    Target<Dim>::pup(p);
    p | name_;
    p | points_;
  }

 private:
  const tnsr::I<DataVector, Dim, Frame::NoFrame>& target_points_no_frame()
      const override {
    return points_;
  }

  std::string name_{"TestTargetOperator"};
  tnsr::I<DataVector, Dim, Frame::NoFrame> points_{};
};

template <size_t Dim>
PUP::able::PUP_ID TestTargetOperator<Dim>::my_PUP_ID = 0;  // NOLINT

template <size_t Dim>
struct TestMetavars {
  using const_global_cache_tags = tmpl::list<domain::Tags::Domain<Dim>>;
  struct factory_creation {
    using factory_classes = tmpl::map<tmpl::pair<
        Target<Dim>, tmpl::list<TestTarget<Dim>, TestTargetOperator<Dim>>>>;
  };

  using component_list = tmpl::list<>;
};

template <size_t Dim>
void test() {
  register_classes_with_charm(
      tmpl::list<TestTarget<Dim>, TestTargetOperator<Dim>>{});
  const size_t number_of_grid_points = 5;
  const double fill = 4.6;

  std::unique_ptr<Target<Dim>> target_base = serialize_and_deserialize(
      std::make_unique<TestTarget<Dim>>(serialize_and_deserialize(
          TestTarget<Dim>{number_of_grid_points, fill})));

  CHECK(target_base->name() == "TestTarget");

  domain::creators::AlignedLattice<Dim> lattice{
      make_array<Dim, std::vector<double>>(std::vector<double>{0.0, 1.0}),
      make_array<Dim, size_t>(0_st),
      {},
      {},
      {},
      {},
      make_array<Dim, bool>(false)};

  const auto box = db::create<tmpl::list<>>();
  Parallel::GlobalCache<TestMetavars<Dim>> cache{{lattice.create_domain()}};

  const auto coords =
      target_base->block_logical_coordinates(box, cache, 1.3, "Grid");

  std::unique_ptr<Target<Dim>> target_operator_base = serialize_and_deserialize(
      std::make_unique<TestTargetOperator<Dim>>(serialize_and_deserialize(
          TestTargetOperator<Dim>{number_of_grid_points, fill})));

  CHECK(target_operator_base->name() == "TestTargetOperator");

  CHECK_THROWS_WITH(
      (target_operator_base->block_logical_coordinates(box, cache, 1.3,
                                                       "Grid")),
      Catch::Matchers::ContainsSubstring("Called the correct operator()"));
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.RuntimeTargets.TargetBase",
    "[Unit]") {
  test<1>();
  test<2>();
  test<3>();
}
}  // namespace intrp::Targets
