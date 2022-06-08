// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Shell.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/ElementInitInterpPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/ElementReceiveInterpPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolationTargetSendPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolationTargetDetail.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/LineSegment.hpp"
#include "Time/Tags.hpp"
#include "Utilities/ProtocolHelpers.hpp"

namespace {

// Simple Variables tag for test.
namespace Tags {
struct TestSolution : db::SimpleTag {
  using type = Scalar<DataVector>;
};
}  // namespace Tags

template <typename Metavariables>
struct mock_element {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = size_t;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename Metavariables::Phase, Metavariables::Phase::Initialization,
      tmpl::list<Actions::SetupDataBox,
                 intrp::Actions::ElementInitInterpPoints<
                     intrp::Tags::InterpPointInfo<Metavariables>>>>>;
  using initial_databox = db::compute_databox_type<
      tmpl::list<intrp::Tags::InterpPointInfo<Metavariables>>>;
  using component_being_mocked =
      DgElementArray<Metavariables, phase_dependent_action_list>;
};

template <typename Metavariables, typename InterpolationTargetTag>
struct mock_interpolation_target {
  static_assert(
      tt::assert_conforms_to<InterpolationTargetTag,
                             intrp::protocols::InterpolationTargetTag>);
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = size_t;
  using const_global_cache_tags = tmpl::push_back<
      Parallel::get_const_global_cache_tags_from_actions<
          tmpl::list<typename InterpolationTargetTag::compute_target_points>>,
      domain::Tags::Domain<Metavariables::volume_dim>>;
  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<Actions::SetupDataBox,
                     intrp::Actions::InitializeInterpolationTarget<
                         Metavariables, InterpolationTargetTag>>>,
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Registration,
          tmpl::list<intrp::Actions::InterpolationTargetSendPointsToElements<
              InterpolationTargetTag>>>>;
  using component_being_mocked =
      intrp::InterpolationTarget<Metavariables, InterpolationTargetTag>;
};

struct MockMetavariables {
  struct InterpolationTargetA
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::TimeStepId;
    using vars_to_interpolate_to_target = tmpl::list<Tags::TestSolution>;
    using compute_items_on_target = tmpl::list<>;
    using compute_target_points =
        ::intrp::TargetPoints::LineSegment<InterpolationTargetA, 3>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<tmpl::list<>,
                                                     InterpolationTargetA>;
    template <typename Metavariables>
    using interpolating_component = mock_element<Metavariables>;
  };
  struct InterpolationTargetB
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::TimeStepId;
    using vars_to_interpolate_to_target = tmpl::list<Tags::TestSolution>;
    using compute_items_on_target = tmpl::list<>;
    using compute_target_points =
        ::intrp::TargetPoints::LineSegment<InterpolationTargetB, 3>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<tmpl::list<>,
                                                     InterpolationTargetA>;
    template <typename Metavariables>
    using interpolating_component = mock_element<Metavariables>;
  };
  static constexpr size_t volume_dim = 3;
  using interpolation_target_tags =
      tmpl::list<InterpolationTargetA, InterpolationTargetB>;

  using component_list = tmpl::list<
      mock_interpolation_target<MockMetavariables, InterpolationTargetA>,
      mock_interpolation_target<MockMetavariables, InterpolationTargetB>,
      mock_element<MockMetavariables>>;
  enum class Phase { Initialization, Registration, Testing, Exit };
};

template <bool IsTimeDep>
void run_test() {
  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();

  using metavars = MockMetavariables;
  using target_component_a =
      mock_interpolation_target<metavars,
                                typename metavars::InterpolationTargetA>;
  using target_component_b =
      mock_interpolation_target<metavars,
                                typename metavars::InterpolationTargetB>;
  using elem_component = mock_element<metavars>;

  // Options
  domain::creators::Shell domain_creator{};
  if constexpr (IsTimeDep) {
    std::unique_ptr<domain::creators::time_dependence::TimeDependence<3>>
        time_dependence = std::make_unique<
            domain::creators::time_dependence::UniformTranslation<3>>(
            0.0, std::array<double, 3>{{0.0, 0.0, 0.0}});
    domain_creator =
        domain::creators::Shell(0.9, 4.9, 1, {{5, 5}}, false, {}, {},
                                {domain::CoordinateMaps::Distribution::Linear},
                                ShellWedges::All, std::move(time_dependence));
  } else {
    domain_creator = domain::creators::Shell(0.9, 4.9, 1, {{5, 5}}, false);
  }

  intrp::OptionHolders::LineSegment<3> line_segment_opts_a(
      {{1.0, 1.0, 1.0}}, {{2.4, 2.4, 2.4}}, 15);
  intrp::OptionHolders::LineSegment<3> line_segment_opts_b(
      {{1.0, 1.0, 1.0}}, {{2.1, 2.1, 2.1}}, 12);
  tuples::TaggedTuple<intrp::Tags::LineSegment<metavars::InterpolationTargetA,
                                               metavars::volume_dim>,
                      domain::Tags::Domain<metavars::volume_dim>,
                      intrp::Tags::LineSegment<metavars::InterpolationTargetB,
                                               metavars::volume_dim>>
      tuple_of_opts{std::move(line_segment_opts_a),
                    domain_creator.create_domain(),
                    std::move(line_segment_opts_b)};

  // Initialization. It's ok if the time independent test has functions of time.
  // They don't get used in either test.
  ActionTesting::MockRuntimeSystem<metavars> runner{
      std::move(tuple_of_opts), {domain_creator.functions_of_time()}};
  ActionTesting::set_phase(make_not_null(&runner),
                           metavars::Phase::Initialization);
  ActionTesting::emplace_component<target_component_a>(&runner, 0);
  for (size_t i = 0; i < 2; ++i) {
    ActionTesting::next_action<target_component_a>(make_not_null(&runner), 0);
  }
  ActionTesting::emplace_component<target_component_b>(&runner, 0);
  for (size_t i = 0; i < 2; ++i) {
    ActionTesting::next_action<target_component_b>(make_not_null(&runner), 0);
  }
  ActionTesting::emplace_component<elem_component>(&runner, 0);
  for (size_t i = 0; i < 2; ++i) {
    ActionTesting::next_action<elem_component>(make_not_null(&runner), 0);
  }
  ActionTesting::set_phase(make_not_null(&runner),
                           metavars::Phase::Registration);

  using point_info_type = std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, metavars::volume_dim,
                                      typename ::Frame::BlockLogical>>>>;

  // element should contain an intrp::Tags::InterpPointInfo<Metavariables>
  // that has an empty point_info for InterpolationTargetA and
  // InterpolationTargetB. This tests ElementInitInterpPoints.
  const auto& init_point_infos =
      ActionTesting::get_databox_tag<elem_component,
                                     ::intrp::Tags::InterpPointInfo<metavars>>(
          runner, 0);

  CHECK(get<intrp::Vars::PointInfoTag<metavars::InterpolationTargetA,
                                      metavars::volume_dim>>(init_point_infos)
            .empty());
  CHECK(get<intrp::Vars::PointInfoTag<metavars::InterpolationTargetB,
                                      metavars::volume_dim>>(init_point_infos)
            .empty());

  // Now invoke the only Registration action (InterpolationTargetSendPoints).
  ActionTesting::next_action<target_component_a>(make_not_null(&runner), 0);
  ActionTesting::next_action<target_component_b>(make_not_null(&runner), 0);
  if constexpr (IsTimeDep) {
    // In this case there should be no simple actions because we didn't send any
    // points to the interpolator. So we just end here because there's nothing
    // else left to do
    CHECK(
        ActionTesting::is_simple_action_queue_empty<elem_component>(runner, 0));
    return;
  }
  // ... and the only queued simple action (ElementReceiveInterpPoints).
  runner.invoke_queued_simple_action<elem_component>(0);
  runner.invoke_queued_simple_action<elem_component>(0);

  ActionTesting::set_phase(make_not_null(&runner), metavars::Phase::Testing);

  // After registration, element should contain an
  // intrp::Tags::InterpPointInfo<Metavariables> that has a
  // point_info for InterpolationTargetA/B containing the
  // expected point info.
  const auto expected_point_info_a = [&domain_creator]() {
    const size_t n_pts = 15;
    tnsr::I<DataVector, 3, Frame::Inertial> points(n_pts);
    for (size_t d = 0; d < 3; ++d) {
      for (size_t i = 0; i < n_pts; ++i) {
        points.get(d)[i] = 1.0 + 0.1 * i;  // Worked out by hand.
      }
    }
    return block_logical_coordinates(domain_creator.create_domain(), points);
  }();
  const auto expected_point_info_b = [&domain_creator]() {
    const size_t n_pts = 12;
    tnsr::I<DataVector, 3, Frame::Inertial> points(n_pts);
    for (size_t d = 0; d < 3; ++d) {
      for (size_t i = 0; i < n_pts; ++i) {
        points.get(d)[i] = 1.0 + 0.1 * i;  // Worked out by hand.
      }
    }
    return block_logical_coordinates(domain_creator.create_domain(), points);
  }();

  const auto& point_infos =
      ActionTesting::get_databox_tag<elem_component,
                                     ::intrp::Tags::InterpPointInfo<metavars>>(
          runner, 0);
  const auto& point_info_a =
      get<intrp::Vars::PointInfoTag<metavars::InterpolationTargetA,
                                    metavars::volume_dim>>(point_infos);
  const auto& point_info_b =
      get<intrp::Vars::PointInfoTag<metavars::InterpolationTargetB,
                                    metavars::volume_dim>>(point_infos);

  auto check_point_info = [](const point_info_type& result,
                             const point_info_type& expected) {
    const size_t number_of_points = expected.size();
    CHECK(result.size() == number_of_points);
    for (size_t i = 0; i < number_of_points; ++i) {
      CHECK(result[i].value().id == expected[i].value().id);
      CHECK_ITERABLE_APPROX(result[i].value().data, expected[i].value().data);
    }
  };

  check_point_info(point_info_a, expected_point_info_a);
  check_point_info(point_info_b, expected_point_info_b);

  // Should be no queued simple actions on either component.
  CHECK(runner.is_simple_action_queue_empty<target_component_a>(0));
  CHECK(runner.is_simple_action_queue_empty<target_component_b>(0));
  CHECK(runner.is_simple_action_queue_empty<elem_component>(0));
}
SPECTRE_TEST_CASE("Unit.NumericalAlgorithms.Interpolator.ElementReceivePoints",
                  "[Unit]") {
  run_test<false>();
  run_test<true>();
}
}  // namespace
