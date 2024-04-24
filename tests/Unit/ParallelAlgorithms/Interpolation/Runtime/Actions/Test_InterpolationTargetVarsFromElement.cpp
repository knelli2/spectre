// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <unordered_set>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/Distribution.hpp"
#include "Domain/Creators/Interval.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Framework/ActionTesting.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/AccessWrapper.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Actions/InterpolationTargetVarsFromElement.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Points/LineSegment.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct UpdateFoT {
  static void apply(
      gsl::not_null<std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>*>
          f_of_t_list,
      const std::string& f_of_t_name, double update_time,
      DataVector update_deriv, double new_expiration_time) {
    (*f_of_t_list)
        .at(f_of_t_name)
        ->update(update_time, std::move(update_deriv), new_expiration_time);
  }
};

namespace Tags {
struct Time : db::SimpleTag {
  using type = double;
};

template <size_t Index>
struct ScalarTag : db::SimpleTag {
  using type = Scalar<DataVector>;
};

template <typename Tag>
struct Sum : db::SimpleTag {
  using type = double;
};

template <typename Tag>
struct SumCompute : Sum<Tag>, db::ComputeTag {
  using base = Sum<Tag>;
  using return_type = double;
  using argument_tags = tmpl::list<Tag>;

  static void function(const gsl::not_null<double*> result,
                       const Scalar<DataVector>& tag) {
    *result = sum(get(tag));
  }
};
}  // namespace Tags

template <typename Target>
struct MockCallback : public intrp2::callbacks::Callback<Target>,
                      public tt::ConformsTo<intrp2::protocols::Callback> {
  using tags_to_observe_on_target =
      tmpl::list<Tags::ScalarTag<0>, Tags::SumCompute<Tags::ScalarTag<0>>>;
  using non_observation_tags_on_target = tmpl::list<>;
  using volume_compute_tags = tmpl::list<>;

  MockCallback() = default;
  explicit MockCallback(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(intrp2::callbacks::Callback<Target>,
                                     MockCallback);
  MockCallback(std::unordered_set<std::string> values_to_observe,
               const size_t expected_num_points)
      : values_to_observe_(std::move(values_to_observe)),
        expected_num_points_(expected_num_points) {}

  using options = tmpl::list<>;
  static constexpr Options::String help = {"halp"};

  void pup(PUP::er& p) override { p | values_to_observe_; }

  const std::unordered_set<std::string>& observables() const override {
    return values_to_observe_;
  }

  template <typename Metavariables>
  void apply(const db::Access& access,
             Parallel::GlobalCache<Metavariables>& cache, const double time) {
    const auto& scalar = db::get<Tags::ScalarTag<0>>(access);
    CHECK(get(scalar).size() == expected_num_points_);
    CHECK(db::get<Tags::Sum<Tags::ScalarTag<0>>>(access) == 0.0);
    (void)cache;
    (void)time;
  }

 private:
  std::unordered_set<std::string> values_to_observe_{};
  size_t expected_num_points_{};
};

template <typename Target>
// NOLINTNEXTLINE
PUP::able::PUP_ID MockCallback<Target>::my_PUP_ID = 0;

using Points = intrp2::points::LineSegment<1>;

struct MockTarget : tt::ConformsTo<intrp2::protocols::Target> {
  using temporal_id_tag = Tags::Time;
  using frame = NoSuchType;
  using points = Points;
  using possible_runtime_callbacks = tmpl::list<MockCallback<MockTarget>>;
  using compile_time_callbacks = tmpl::list<>;
};

template <typename Metavariables>
struct MockTargetComponent {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = std::string;
  using const_global_cache_tags = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 tmpl::list<intrp2::Tags::AccessWrapper>>>>>;
};

struct TestMetavars {
  static constexpr size_t volume_dim = 1;

  using component_list = tmpl::list<MockTargetComponent<TestMetavars>>;

  using const_global_cache_tags = tmpl::list<domain::Tags::Domain<volume_dim>>;
  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize>;

  struct factory_creation {
    using factory_classes = tmpl::map<
        tmpl::pair<intrp2::callbacks::Callback<MockTarget>,
                   tmpl::list<MockCallback<MockTarget>>>,
        tmpl::pair<
            intrp2::AccessWrapper,
            tmpl::list<intrp2::TargetAccessWrapper<MockTarget, volume_dim>>>>;
  };
};

void test() {
  using target_component = MockTargetComponent<TestMetavars>;
  using action = intrp2::Actions::ReceiveVolumeTensors<MockTarget>;
  constexpr size_t Dim = TestMetavars::volume_dim;
  // Registration hell... :(
  register_factory_classes_with_charm<TestMetavars>();
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();

  // This makes the points nicely spaced every 0.1
  const size_t num_points = 11;
  Points line_segment{{0.0}, {1.0}, num_points};
  std::vector<size_t> all_offsets(num_points);
  std::iota(all_offsets.begin(), all_offsets.end(), 0_st);
  const double invalid_fill = std::numeric_limits<double>::max();
  std::unordered_set<std::string> names_of_vars{"Scalar0", "Sum"};

  std::unique_ptr<intrp2::AccessWrapper> initial_access_wrapper =
      std::make_unique<intrp2::TargetAccessWrapper<MockTarget, Dim>>();

  // This is normally done by the Event, so we have to do it manually here
  {
    db::Access& access = initial_access_wrapper->access();
    db::mutate<intrp2::Tags::Frame, intrp2::Tags::Points<Dim>,
               intrp2::Tags::NumberOfSetsOfPoints,
               intrp2::Tags::InvalidPointsFillValue,
               intrp2::Tags::Callbacks<MockTarget>>(
        [&line_segment, &invalid_fill, &names_of_vars, &num_points](
            const gsl::not_null<std::string*> frame,
            const gsl::not_null<tnsr::I<DataVector, Dim, Frame::NoFrame>*>
                points,
            const gsl::not_null<size_t*> number_of_sets_of_points,
            const gsl::not_null<double*> invalid_fill_value,
            const gsl::not_null<std::vector<
                std::unique_ptr<intrp2::callbacks::Callback<MockTarget>>>*>
                callbacks) {
          *frame = "Inertial";
          *points = line_segment.target_points_no_frame();
          *number_of_sets_of_points = 1;
          *invalid_fill_value = invalid_fill;
          callbacks->resize(1);
          (*callbacks)[0] = std::make_unique<MockCallback<MockTarget>>(
              names_of_vars, num_points);
        },
        make_not_null(&access));
  }
  // This is a compute tag, so we only send over its argument tags
  names_of_vars.erase("Sum");

  const double temporal_id_1 = 1.0;
  const double temporal_id_2 = 1.3;
  const std::string array_index = pretty_type::name<MockTarget>();
  // We give the domain a velocity so that at temporal_id_1 all points are
  // valid, but at temporal_id_2 some will be invalid and we can test that they
  // were filled properly.
  const domain::creators::Interval domain_creator{
      {0.0},
      {1.0},
      {3_st},
      {4_st},
      {false},
      domain::CoordinateMaps::Distribution::Linear,
      std::nullopt,
      std::make_unique<
          domain::creators::time_dependence::UniformTranslation<Dim>>(
          temporal_id_1, std::array{1.0})};
  // Even though the time dependence adds a FoT, it never expires. We need it to
  // expire so we can test the checking of the FoTs, so we make a new FoT that
  // does expire between the two times and overwrite the existing one
  auto functions_of_time = domain_creator.functions_of_time();
  const double first_expiration_time = 0.5 * (temporal_id_1 + temporal_id_2);
  {
    std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>
        trans_with_expiration =
            std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<1>>(
                temporal_id_1, std::array{DataVector{0.0}, DataVector{1.0}},
                first_expiration_time);
    functions_of_time["Translation"] = std::move(trans_with_expiration);
  }

  ActionTesting::MockRuntimeSystem<TestMetavars> runner{
      {domain_creator.create_domain()}, {std::move(functions_of_time)}};
  ActionTesting::emplace_array_component_and_initialize<target_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, array_index,
      {std::move(initial_access_wrapper)});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto& cache = ActionTesting::cache<target_component>(runner, array_index);
  auto& box = ActionTesting::get_databox<target_component>(
      make_not_null(&runner), array_index);
  auto& access_wrapper = db::get_mutable_reference<intrp2::Tags::AccessWrapper>(
      make_not_null(&box));
  auto& access = access_wrapper->access();

  bool check_for_invalid_points = false;
  std::unordered_map<std::string, std::vector<DataVector>> interpolated_vars{};
  std::vector<size_t> global_offsets{0, 8, 9};
  std::unordered_set<size_t> cumulative_global_offsets_1{};
  std::unordered_set<size_t> cumulative_global_offsets_2{};

  const auto fill_interpolated_vars =
      [&names_of_vars, &interpolated_vars,
       &global_offsets](std::unordered_set<size_t>& cumulative_global_offsets) {
        cumulative_global_offsets.insert(global_offsets.begin(),
                                         global_offsets.end());
        for (const std::string& tensor_name : names_of_vars) {
          // For now each tensor is a scalar so there's only cone component
          interpolated_vars[tensor_name] = std::vector<DataVector>(1);
          // TODO: Change from zeros to better data
          interpolated_vars.at(tensor_name)[0] =
              DataVector{global_offsets.size(), 0.0};
        }
      };

  fill_interpolated_vars(cumulative_global_offsets_1);

  // Nothing should be initialized yet
  {
    INFO("Checking that nothing is initialized yet");
    const auto check = [&access](const double temporal_id) {
      CHECK_FALSE(
          db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
              temporal_id));
      CHECK_FALSE(
          db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access).contains(
              temporal_id));
      CHECK_FALSE(
          db::get<intrp2::Tags::InterpolatedVars<double>>(access).contains(
              temporal_id));
      CHECK_FALSE(
          db::get<intrp2::Tags::InvalidIndices<double>>(access).contains(
              temporal_id));
      CHECK(
          db::get<intrp2::Tags::CompletedTemporalIds<double>>(access).empty());
    };
    check(temporal_id_1);
    check(temporal_id_2);
  }

  // Send the data to the target at temporal_id_1
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_1, interpolated_vars,
      global_offsets, check_for_invalid_points);

  // Check that everything was received correctly at temporal_id_1
  {
    INFO("Checking first data at temporal_id_1");
    CHECK(db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
        temporal_id_1));
    const auto& number_of_filled_points =
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access);
    CHECK(number_of_filled_points.contains(temporal_id_1));
    CHECK(number_of_filled_points.at(temporal_id_1) == 3);
    const auto& stored_vars =
        db::get<intrp2::Tags::InterpolatedVars<double>>(access);
    CHECK(stored_vars.contains(temporal_id_1));
    const auto& stored_vars_at_temporal_id_1 = stored_vars.at(temporal_id_1);
    for (const std::string& tensor_name : names_of_vars) {
      CHECK(stored_vars_at_temporal_id_1.contains(tensor_name));
      CHECK(stored_vars_at_temporal_id_1.at(tensor_name).size() == 1);
      for (size_t global_offset : all_offsets) {
        // TODO: Needs to be updated when we switch from zeros to real data
        CHECK(stored_vars_at_temporal_id_1.at(tensor_name)[0][global_offset] ==
              (alg::found(cumulative_global_offsets_1, global_offset)
                   ? 0.0
                   : invalid_fill));
      }
    }
    CHECK_FALSE(db::get<intrp2::Tags::InvalidIndices<double>>(access).contains(
        temporal_id_1));
    CHECK(db::get<intrp2::Tags::CompletedTemporalIds<double>>(access).empty());
  }

  interpolated_vars.clear();
  global_offsets = std::vector<size_t>{3, 4, 5};
  fill_interpolated_vars(cumulative_global_offsets_2);

  // Send the data to the target at temporal_id_2. Also check for invalid
  // points.
  check_for_invalid_points = true;
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_2, interpolated_vars,
      global_offsets, check_for_invalid_points);
  check_for_invalid_points = false;

  // temporal_id_2 is after the expiration of the FoT so we should have moved
  // the data into the target component, but we shouldn't have checked the
  // invalid points yet.
  {
    INFO("Checking first data at temporal_id_2 with FoT callback queued");
    CHECK(db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
        temporal_id_2));
    const auto& number_of_filled_points =
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access);
    CHECK(number_of_filled_points.contains(temporal_id_2));
    CHECK(number_of_filled_points.at(temporal_id_2) == 3);
    const auto& stored_vars =
        db::get<intrp2::Tags::InterpolatedVars<double>>(access);
    CHECK(stored_vars.contains(temporal_id_2));
    const auto& stored_vars_at_temporal_id_2 = stored_vars.at(temporal_id_2);
    for (const std::string& tensor_name : names_of_vars) {
      CHECK(stored_vars_at_temporal_id_2.contains(tensor_name));
      CHECK(stored_vars_at_temporal_id_2.at(tensor_name).size() == 1);
      for (size_t global_offset : all_offsets) {
        // TODO: Needs to be updated when we switch from zeros to real data
        CHECK(stored_vars_at_temporal_id_2.at(tensor_name)[0][global_offset] ==
              (alg::found(cumulative_global_offsets_2, global_offset)
                   ? 0.0
                   : invalid_fill));
      }
    }
    const auto& invalid_points =
        db::get<intrp2::Tags::InvalidIndices<double>>(access);
    CHECK_FALSE(invalid_points.contains(temporal_id_2));
  }

  // Now update the FoT
  Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
      cache, "Translation"s, first_expiration_time, DataVector{1.0},
      temporal_id_2 + 1.0);

  // There should be a simple action callback queued. Check that there is and
  // call it.
  CHECK(ActionTesting::number_of_queued_simple_actions<target_component>(
            runner, array_index) == 1);
  ActionTesting::invoke_queued_simple_action<target_component>(
      make_not_null(&runner), array_index);

  // Now we will actually check if there are any invalid points
  {
    INFO("Checking invalid points at temporal_id_2");
    const auto& invalid_points =
        db::get<intrp2::Tags::InvalidIndices<double>>(access);
    CHECK(invalid_points.contains(temporal_id_2));
    // These are the indices which should be outside the domain now
    CHECK(invalid_points.at(temporal_id_2) ==
          std::unordered_set<size_t>{0, 1, 2});
    CHECK(db::get<intrp2::Tags::CompletedTemporalIds<double>>(access).empty());
  }

  interpolated_vars.clear();
  global_offsets = std::vector<size_t>{1, 5, 6, 7};
  fill_interpolated_vars(cumulative_global_offsets_1);

  // Second round of data for temporal_id_1. This time check for invalid points
  // as well
  check_for_invalid_points = true;
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_1, interpolated_vars,
      global_offsets, check_for_invalid_points);
  check_for_invalid_points = false;

  // Check that everything was received correctly
  {
    INFO("Checking second data and invalid points at temporal_id_1.");
    CHECK(db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
        temporal_id_1));
    const auto& number_of_filled_points =
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access);
    CHECK(number_of_filled_points.contains(temporal_id_1));
    CHECK(number_of_filled_points.at(temporal_id_1) == 7);
    const auto& stored_vars =
        db::get<intrp2::Tags::InterpolatedVars<double>>(access);
    CHECK(stored_vars.contains(temporal_id_1));
    const auto& stored_vars_at_temporal_id_1 = stored_vars.at(temporal_id_1);
    for (const std::string& tensor_name : names_of_vars) {
      CHECK(stored_vars_at_temporal_id_1.contains(tensor_name));
      CHECK(stored_vars_at_temporal_id_1.at(tensor_name).size() == 1);
      for (size_t global_offset : all_offsets) {
        // TODO: Needs to be updated when we switch from zeros to real data
        CHECK(stored_vars_at_temporal_id_1.at(tensor_name)[0][global_offset] ==
              (alg::found(cumulative_global_offsets_1, global_offset)
                   ? 0.0
                   : invalid_fill));
      }
    }
    const auto& invalid_points =
        db::get<intrp2::Tags::InvalidIndices<double>>(access);
    CHECK(invalid_points.contains(temporal_id_1));
    // We don't actually have invalid points for this temporal id
    CHECK(invalid_points.at(temporal_id_1).empty());
    CHECK(db::get<intrp2::Tags::CompletedTemporalIds<double>>(access).empty());
  }

  interpolated_vars.clear();
  global_offsets = std::vector<size_t>{8, 9, 10};
  fill_interpolated_vars(cumulative_global_offsets_2);

  // Second round of data for temporal_id_2.
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_2, interpolated_vars,
      global_offsets, check_for_invalid_points);

  {
    INFO("Checking second data at temporal_id_2");
    CHECK(db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
        temporal_id_2));
    const auto& number_of_filled_points =
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access);
    CHECK(number_of_filled_points.contains(temporal_id_2));
    CHECK(number_of_filled_points.at(temporal_id_2) == 6);
    const auto& stored_vars =
        db::get<intrp2::Tags::InterpolatedVars<double>>(access);
    CHECK(stored_vars.contains(temporal_id_2));
    const auto& stored_vars_at_temporal_id_2 = stored_vars.at(temporal_id_2);
    for (const std::string& tensor_name : names_of_vars) {
      CHECK(stored_vars_at_temporal_id_2.contains(tensor_name));
      CHECK(stored_vars_at_temporal_id_2.at(tensor_name).size() == 1);
      for (size_t global_offset : all_offsets) {
        // TODO: Needs to be updated when we switch from zeros to real data
        CHECK(stored_vars_at_temporal_id_2.at(tensor_name)[0][global_offset] ==
              (alg::found(cumulative_global_offsets_2, global_offset)
                   ? 0.0
                   : invalid_fill));
      }
    }
    const auto& invalid_points =
        db::get<intrp2::Tags::InvalidIndices<double>>(access);
    CHECK(invalid_points.contains(temporal_id_2));
    CHECK(invalid_points.at(temporal_id_2) ==
          std::unordered_set<size_t>{0, 1, 2});
    CHECK(db::get<intrp2::Tags::CompletedTemporalIds<double>>(access).empty());
  }

  interpolated_vars.clear();
  global_offsets = std::vector<size_t>{2, 3, 4, 10};
  fill_interpolated_vars(cumulative_global_offsets_1);

  // Final round of data for temporal_id_1. This should trigger the callbacks
  // and clean up everything
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_1, interpolated_vars,
      global_offsets, check_for_invalid_points);

  {
    INFO("Checking final data at temporal_id_1");
    CHECK_FALSE(
        db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
            temporal_id_1));
    CHECK_FALSE(
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access).contains(
            temporal_id_1));
    CHECK_FALSE(
        db::get<intrp2::Tags::InterpolatedVars<double>>(access).contains(
            temporal_id_1));
    CHECK_FALSE(db::get<intrp2::Tags::InvalidIndices<double>>(access).contains(
        temporal_id_1));
    const auto& completed_temporal_ids =
        db::get<intrp2::Tags::CompletedTemporalIds<double>>(access);
    CHECK(completed_temporal_ids.size() == 1);
    CHECK(completed_temporal_ids.front() == temporal_id_1);
  }

  interpolated_vars.clear();
  global_offsets = std::vector<size_t>{6, 7};
  fill_interpolated_vars(cumulative_global_offsets_2);

  // Final round of data for temporal_id_2. This should trigger the callbacks
  // and clean up everything. Note that we didn't send all the points because
  // some were invalid
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_2, interpolated_vars,
      global_offsets, check_for_invalid_points);

  {
    INFO("Checking final data at temporal_id_2");
    CHECK_FALSE(
        db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
            temporal_id_2));
    CHECK_FALSE(
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access).contains(
            temporal_id_2));
    CHECK_FALSE(
        db::get<intrp2::Tags::InterpolatedVars<double>>(access).contains(
            temporal_id_2));
    CHECK_FALSE(db::get<intrp2::Tags::InvalidIndices<double>>(access).contains(
        temporal_id_2));
    const auto& completed_temporal_ids =
        db::get<intrp2::Tags::CompletedTemporalIds<double>>(access);
    CHECK(completed_temporal_ids.size() == 2);
    CHECK(completed_temporal_ids.front() == temporal_id_2);
  }

  // Check that sending more data at the same temporal_id_1 does nothing
  ActionTesting::simple_action<target_component, action>(
      make_not_null(&runner), array_index, temporal_id_1, interpolated_vars,
      global_offsets, check_for_invalid_points);

  {
    INFO("Checking extra send for temporal_id_1");
    CHECK_FALSE(
        db::get<intrp2::Tags::CurrentTemporalIds<double>>(access).contains(
            temporal_id_1));
    CHECK_FALSE(
        db::get<intrp2::Tags::NumberOfFilledPoints<double>>(access).contains(
            temporal_id_1));
    CHECK_FALSE(
        db::get<intrp2::Tags::InterpolatedVars<double>>(access).contains(
            temporal_id_1));
    CHECK_FALSE(db::get<intrp2::Tags::InvalidIndices<double>>(access).contains(
        temporal_id_1));
    const auto& completed_temporal_ids =
        db::get<intrp2::Tags::CompletedTemporalIds<double>>(access);
    CHECK(completed_temporal_ids.size() == 2);
    CHECK(completed_temporal_ids.front() == temporal_id_2);
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.Interpolation.Runtime.VarsFromElement", "[Unit]") {
  test();
}
