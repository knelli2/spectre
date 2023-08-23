// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/Tags/SystemTags.hpp"
#include "ControlSystem/UpdateFunctionOfTime.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/SettleToConstant.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Utilities/CloneUniquePtrs.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct System1 {
  static std::string name() { return "FoT1"; }
};
struct System2 {
  static std::string name() { return "FoT2"; }
};
struct System3 {
  static std::string name() { return "FoT3"; }
};

template <typename Metavariables, size_t Index = 1>
struct TestSingleton {
  static std::string name() { return "FoT" + get_output(Index); }
  using chare_type = ActionTesting::MockSingletonChare;
  using array_index = size_t;
  using metavariables = Metavariables;
  using const_global_cache_tags =
      tmpl::list<control_system::Tags::SystemToCombinedNames>;
  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize,
                 control_system::Tags::MeasurementTimescales>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<
          tmpl::list<control_system::Tags::UpdateAggregators>>>>>;
};

struct TestingMetavariables {
  using component_list = tmpl::list<TestSingleton<TestingMetavariables>>;
};

void test_mutates() {
  constexpr size_t deriv_order = 2;
  const double t0 = 0.0;
  const double expr_time = 1.0;
  double update_time = 1.5;
  double expiration_time = 2.0;
  domain::FunctionsOfTime::register_derived_with_charm();

  // Construct unordered map
  const std::string pp_name{"TestPiecewisePolynomial"};
  const std::string quatfot_name{"TestQuatFunctionOfTime"};
  std::array<DataVector, deriv_order + 1> init_pp_func{
      {DataVector{3, 0.0}, DataVector{3, 0.0}, DataVector{3, 0.0}}};
  std::array<DataVector, 1> init_quat_func{DataVector{4, 1.0}};

  domain::FunctionsOfTime::PiecewisePolynomial<deriv_order> expected_pp{
      t0, init_pp_func, expr_time};
  domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>
      expected_quatfot{t0, init_quat_func, init_pp_func, expr_time};

  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      measurement_timescales{};
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      f_of_t_map{};
  f_of_t_map[pp_name] = std::make_unique<
      domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>>(
      t0, init_pp_func, expr_time);
  f_of_t_map[quatfot_name] = std::make_unique<
      domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>>(
      t0, init_quat_func, init_pp_func, expr_time);

  const DataVector updated_deriv{3, 1.0};

  // Construct mock system and set to Testing phase
  ActionTesting::MockRuntimeSystem<TestingMetavariables> runsys{
      {}, {std::move(f_of_t_map), std::move(measurement_timescales)}};
  ActionTesting::emplace_singleton_component_and_initialize<
      TestSingleton<TestingMetavariables>>(make_not_null(&runsys),
                                           ActionTesting::NodeId{0},
                                           ActionTesting::LocalCoreId{0}, {});
  ActionTesting::set_phase(make_not_null(&runsys), Parallel::Phase::Testing);

  // Update functions of time in global cache with new deriv
  auto& cache =
      ActionTesting::cache<TestSingleton<TestingMetavariables>>(runsys, 0_st);

  for (auto& name : {pp_name, quatfot_name}) {
    Parallel::mutate<domain::Tags::FunctionsOfTime,
                     control_system::UpdateSingleFunctionOfTime>(
        cache, name, update_time, updated_deriv, expiration_time);
  }

  // Update expected function of time
  expected_pp.update(update_time, updated_deriv, expiration_time);
  expected_quatfot.update(update_time, updated_deriv, expiration_time);

  // Check that the FunctionsOfTime are what we expected
  const auto& cache_f_of_t_map = get<domain::Tags::FunctionsOfTime>(cache);
  {
    const auto& cache_pp = dynamic_cast<
        const domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>&>(
        *(cache_f_of_t_map.at(pp_name)));
    const auto& cache_quatfot = dynamic_cast<
        const domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>&>(
        *(cache_f_of_t_map.at(quatfot_name)));

    CHECK(cache_pp == expected_pp);
    CHECK(cache_quatfot == expected_quatfot);
  }
  // Update functions of time in global cache with new expiration time
  expiration_time += 1.0;
  for (auto& name : {pp_name, quatfot_name}) {
    Parallel::mutate<domain::Tags::FunctionsOfTime,
                     control_system::ResetFunctionOfTimeExpirationTime>(
        cache, name, expiration_time);
  }

  // Update expected function of time
  expected_pp.reset_expiration_time(expiration_time);
  expected_quatfot.reset_expiration_time(expiration_time);
  {
    const auto& cache_pp = dynamic_cast<
        const domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>&>(
        *(cache_f_of_t_map.at(pp_name)));
    const auto& cache_quatfot = dynamic_cast<
        const domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>&>(
        *(cache_f_of_t_map.at(quatfot_name)));

    // Check that the FunctionsOfTime and expiration times are what we expected
    CHECK(cache_pp == expected_pp);
    CHECK(cache_quatfot == expected_quatfot);
    CHECK(cache_pp.time_bounds()[1] == expiration_time);
    CHECK(cache_quatfot.time_bounds()[1] == expiration_time);
  }

  // Check updating both at the same time
  update_time = expiration_time;
  expiration_time += 1.0;
  std::unordered_map<std::string, std::pair<DataVector, double>> update_args{};
  for (const auto& name : {pp_name, quatfot_name}) {
    update_args[name] = std::make_pair(updated_deriv, expiration_time);
  }

  Parallel::mutate<domain::Tags::FunctionsOfTime,
                   control_system::UpdateMultipleFunctionsOfTime>(
      cache, update_time, std::move(update_args));

  expected_pp.update(update_time, updated_deriv, expiration_time);
  expected_quatfot.update(update_time, updated_deriv, expiration_time);

  {
    const auto& cache_pp = dynamic_cast<
        const domain::FunctionsOfTime::PiecewisePolynomial<deriv_order>&>(
        *(cache_f_of_t_map.at(pp_name)));
    const auto& cache_quatfot = dynamic_cast<
        const domain::FunctionsOfTime::QuaternionFunctionOfTime<deriv_order>&>(
        *(cache_f_of_t_map.at(quatfot_name)));

    CHECK(cache_pp == expected_pp);
    CHECK(cache_quatfot == expected_quatfot);
  }
}

void test_update_aggregator() {
  std::unordered_set<std::string> names{"Illinois", "California", "Wisconsin",
                                        "Florida", "Texas"};
  std::vector<std::string> sorted_names(names.begin(), names.end());
  alg::sort(sorted_names);

  control_system::UpdateAggregator unserialized_aggregator{names};
  control_system::UpdateAggregator aggregator =
      serialize_and_deserialize(unserialized_aggregator);

  CHECK(aggregator.combined_name() ==
        "CaliforniaFloridaIllinoisTexasWisconsin");

  aggregator.insert("Wisconsin", DataVector{1.0}, 1.0, DataVector{5.0}, 5.0);
  CHECK_FALSE(aggregator.is_ready());
  aggregator.insert("Texas", DataVector{2.0}, 2.0, DataVector{4.0}, 4.0);
  CHECK_FALSE(aggregator.is_ready());
  aggregator.insert("Illinois", DataVector{3.0}, 3.0, DataVector{3.0}, 3.0);
  CHECK_FALSE(aggregator.is_ready());
  aggregator.insert("Florida", DataVector{4.0}, 4.0, DataVector{2.0}, 2.0);
  CHECK_FALSE(aggregator.is_ready());
  aggregator.insert("California", DataVector{5.0}, 5.0, DataVector{1.0}, 1.0);
  CHECK(aggregator.is_ready());

  const std::unordered_map<std::string, std::pair<DataVector, double>>
      combined_fot = aggregator.combined_fot_expiration_times();
  CHECK(aggregator.is_ready());
  const std::pair<double, double> combined_measurements =
      aggregator.combined_measurement_expiration_time();
  CHECK_FALSE(aggregator.is_ready());

  CHECK(combined_fot.size() == sorted_names.size());
  for (size_t i = 0; i < sorted_names.size(); i++) {
    CHECK(combined_fot.count(sorted_names[i]) == 1);
    CHECK(combined_fot.at(sorted_names[i]) ==
          std::make_pair(DataVector{double(i + 1)}, 1.0));
  }

  CHECK(combined_measurements == std::make_pair(1.0, 1.0));
}

struct AggregatorMetavariables {
  using metavariables = AggregatorMetavariables;
  using component_list = tmpl::list<TestSingleton<metavariables, 1>,
                                    TestSingleton<metavariables, 2>,
                                    TestSingleton<metavariables, 3>>;
};

void test_aggregate_update_action() {
  domain::FunctionsOfTime::register_derived_with_charm();
  using metavars = AggregatorMetavariables;
  using component1 = TestSingleton<metavars, 1>;
  using component2 = TestSingleton<metavars, 2>;
  using component3 = TestSingleton<metavars, 3>;
  const double t0 = 0.0;
  const double old_fot_expiration13 = 4.0;
  const double old_fot_expiration2 = 2.0;
  const double old_measurement_expiration13 = old_fot_expiration13 - 0.5;
  const double old_measurement_expiration2 = old_fot_expiration2 - 0.5;

  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      f_of_t_map{};
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      measurement_map{};

  // We must set up the functions of time and measurement timescales in a
  // consistent state. Here FoT1 and FoT3 use the same measurement timescale,
  // and this have the same expiration time. FoT2 uses a different measurement
  // timescale so it has it's own expiration
  f_of_t_map["FoT1"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          t0, std::array{DataVector{1.0}}, old_fot_expiration13);
  f_of_t_map["FoT2"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          t0, std::array{DataVector{1.0}}, old_fot_expiration2);
  f_of_t_map["FoT3"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          t0, std::array{DataVector{1.0}}, old_fot_expiration13);

  measurement_map["FoT1FoT3"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          t0, std::array{DataVector{4.0}}, old_measurement_expiration13);
  measurement_map["FoT2"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
          t0, std::array{DataVector{2.0}}, old_measurement_expiration2);

  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      expected_f_of_t_map = clone_unique_ptrs(f_of_t_map);
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      expected_measurement_map = clone_unique_ptrs(measurement_map);

  control_system::UpdateAggregator aggregator13{{"FoT1", "FoT3"}};
  control_system::UpdateAggregator aggregator2{{"FoT2"}};
  std::unordered_map<std::string, control_system::UpdateAggregator>
      aggregators{};
  aggregators["FoT1FoT3"] = aggregator13;
  aggregators["FoT2"] = aggregator2;

  std::unordered_map<std::string, std::string> system_to_combined_names{};
  system_to_combined_names["FoT1"] = "FoT1FoT3";
  system_to_combined_names["FoT2"] = "FoT2";
  system_to_combined_names["FoT3"] = "FoT1FoT3";

  ActionTesting::MockRuntimeSystem<metavars> runner{
      {std::move(system_to_combined_names)},
      {std::move(f_of_t_map), std::move(measurement_map)}};
  ActionTesting::emplace_singleton_component_and_initialize<component1>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, {std::move(aggregators)});
  ActionTesting::emplace_singleton_component_and_initialize<component2>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, {});
  ActionTesting::emplace_singleton_component_and_initialize<component3>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, {});
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  const auto& cache = ActionTesting::cache<component1>(runner, 0_st);
  const auto& functions_of_time =
      Parallel::get<domain::Tags::FunctionsOfTime>(cache);
  const auto& measurement_timescales =
      Parallel::get<control_system::Tags::MeasurementTimescales>(cache);

  using tag = control_system::Tags::UpdateAggregators;
  const std::unordered_map<std::string, control_system::UpdateAggregator>&
      box_aggregators_component1 =
          ActionTesting::get_databox_tag<component1, tag>(runner, 0);
  const std::unordered_map<std::string, control_system::UpdateAggregator>&
      box_aggregators_component2 =
          ActionTesting::get_databox_tag<component2, tag>(runner, 0);
  const std::unordered_map<std::string, control_system::UpdateAggregator>&
      box_aggregators_component3 =
          ActionTesting::get_databox_tag<component3, tag>(runner, 0);

  CHECK_FALSE(box_aggregators_component1.empty());
  CHECK(box_aggregators_component2.empty());
  CHECK(box_aggregators_component3.empty());

  const control_system::UpdateAggregator& box_aggregator13 =
      box_aggregators_component1.at("FoT1FoT3");
  const control_system::UpdateAggregator& box_aggregator2 =
      box_aggregators_component1.at("FoT2");

  //   const double time = 1.0;
  const DataVector new_measurement_timescale1{2.0};
  const DataVector new_measurement_timescale2{3.0};
  const DataVector new_measurement_timescale3{4.0};
  const DataVector new_measurement_timescale13{
      std::min(new_measurement_timescale1[0], new_measurement_timescale3[0])};
  const double new_fot_expiration1 = old_fot_expiration13 + 2.0;
  const double new_fot_expiration2 = old_fot_expiration2 + 1.0;
  const double new_fot_expiration3 = old_fot_expiration13 + 1.5;
  const double new_fot_expiration13 =
      std::min(new_fot_expiration1, new_fot_expiration3);
  const DataVector control_signal{1.0};
  const double new_measurement_expiration1 = new_fot_expiration1 - 0.5;
  const double new_measurement_expiration2 = new_fot_expiration2 - 0.5;
  const double new_measurement_expiration3 = new_fot_expiration3 - 0.5;
  const double new_measurement_expiration13 =
      std::min(new_measurement_expiration1, new_measurement_expiration3);

  // The aggregators will never "be ready" to use because when they actually are
  // ready, they are reset immediately during the simple action
  CHECK_FALSE(box_aggregator13.is_ready());
  CHECK_FALSE(box_aggregator2.is_ready());
  CHECK(box_aggregators_component2.empty());
  CHECK(box_aggregators_component3.empty());

  // Update one of the two functions of time for the 13 measurement
  ActionTesting::simple_action<component1,
                               control_system::AggregateUpdate<System1>>(
      make_not_null(&runner), /*time,*/ 0_st, new_measurement_timescale1,
      old_measurement_expiration13, new_measurement_expiration1, control_signal,
      old_fot_expiration13, new_fot_expiration1);

  CHECK_FALSE(box_aggregator13.is_ready());
  CHECK_FALSE(box_aggregator2.is_ready());
  CHECK(box_aggregators_component2.empty());
  CHECK(box_aggregators_component3.empty());

  // Update the 2 measurement and check that the function of time and
  // measurement timescale have been updated appropriately
  ActionTesting::simple_action<component1,
                               control_system::AggregateUpdate<System2>>(
      make_not_null(&runner), /*time,*/ 0_st, new_measurement_timescale2,
      old_measurement_expiration2, new_measurement_expiration2, control_signal,
      old_fot_expiration2, new_fot_expiration2);

  CHECK_FALSE(box_aggregator13.is_ready());
  CHECK_FALSE(box_aggregator2.is_ready());
  CHECK(box_aggregators_component2.empty());
  CHECK(box_aggregators_component3.empty());

  expected_f_of_t_map["FoT2"]->update(old_fot_expiration2, control_signal,
                                      new_fot_expiration2);
  expected_measurement_map["FoT2"]->update(old_measurement_expiration2,
                                           new_measurement_timescale2,
                                           new_measurement_expiration2);

  const auto check_equal =
      [](const std::string& name,
         const std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>& test,
         const std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
             expected) {
        CAPTURE(name);
        CHECK(dynamic_cast<
                  const domain::FunctionsOfTime::PiecewisePolynomial<0>&>(
                  *test.at(name)) ==
              dynamic_cast<
                  const domain::FunctionsOfTime::PiecewisePolynomial<0>&>(
                  *expected.at(name)));
      };

  check_equal("FoT2", functions_of_time, expected_f_of_t_map);
  check_equal("FoT2", measurement_timescales, expected_measurement_map);

  // Update the second of the 13 measurement which should update both functions
  // of time and just the one measurement timescale
  ActionTesting::simple_action<component1,
                               control_system::AggregateUpdate<System3>>(
      make_not_null(&runner), /*time,*/ 0_st, new_measurement_timescale3,
      old_measurement_expiration13, new_measurement_expiration3, control_signal,
      old_fot_expiration13, new_fot_expiration3);

  CHECK_FALSE(box_aggregator13.is_ready());
  CHECK_FALSE(box_aggregator2.is_ready());
  CHECK(box_aggregators_component2.empty());
  CHECK(box_aggregators_component3.empty());

  expected_f_of_t_map["FoT1"]->update(old_fot_expiration13, control_signal,
                                      new_fot_expiration13);
  expected_f_of_t_map["FoT3"]->update(old_fot_expiration13, control_signal,
                                      new_fot_expiration13);
  expected_measurement_map["FoT1FoT3"]->update(old_measurement_expiration13,
                                               new_measurement_timescale13,
                                               new_measurement_expiration13);

  check_equal("FoT1", functions_of_time, expected_f_of_t_map);
  check_equal("FoT3", functions_of_time, expected_f_of_t_map);
  check_equal("FoT1FoT3", measurement_timescales, expected_measurement_map);
}

SPECTRE_TEST_CASE("Unit.ControlSystem.UpdateFunctionOfTime",
                  "[Unit][ControlSystem]") {
  test_mutates();
  test_update_aggregator();
  test_aggregate_update_action();
}
}  // namespace
