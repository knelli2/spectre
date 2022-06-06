// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <deque>
#include <functional>
#include <memory>
#include <pup.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/IdPair.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Block.hpp"
#include "Domain/BlockLogicalCoordinates.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Shell.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/LogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolationTarget.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorReceiveVolumeData.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/Interpolation/Actions/InterpolatorRegisterElement.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/Interpolation/Actions/TryToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/ObserveTimeSeriesOnSurface.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolatedVars.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/ComputeVarsToInterpolate.hpp"
#include "ParallelAlgorithms/Interpolation/Protocols/InterpolationTargetTag.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/LineSegment.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Rational.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma: no_forward_declare ActionTesting::InitializeDataBox
// IWYU pragma: no_forward_declare Tensor
namespace intrp::Actions {
template <typename InterpolationTargetTag>
struct AddTemporalIdsToInterpolationTarget;
template <typename InterpolationTargetTag>
struct InterpolationTargetReceiveVars;
}  // namespace intrp::Actions
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace db {
template <typename TagsList>
class DataBox;
}  // namespace db
namespace intrp::Tags {
template <typename TemporalId>
struct InterpolatedVarsHolders;
template <typename TemporalId>
struct TemporalIds;
template <typename Metavariables, typename TemporalId>
struct VolumeVarsInfo;
}  // namespace intrp::Tags

namespace {

// Simple DataBoxItems for test.
namespace Tags {
struct Square : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct SquareCompute : Square, db::ComputeTag {
  static void function(gsl::not_null<Scalar<DataVector>*> result,
                       const Scalar<DataVector>& x) {
    get(*result) = square(get(x));
  }
  using argument_tags = tmpl::list<gr::Tags::Lapse<DataVector>>;
  using base = Square;
  using return_type = Scalar<DataVector>;
};
}  // namespace Tags

struct ComputeSquare
    : tt::ConformsTo<intrp::protocols::ComputeVarsToInterpolate> {
  template <typename SrcTag, typename DestTag>
  static void apply(
      const gsl::not_null<Variables<tmpl::list<DestTag>>*> target_vars,
      const Variables<tmpl::list<SrcTag>>& src_vars,
      const Mesh<3>& /* mesh */) {
    get(get<DestTag>(*target_vars)) = square(get(get<SrcTag>(src_vars)));
  }

  using allowed_src_tags = tmpl::list<>;
  using required_src_tags = tmpl::list<>;
  template <typename Frame>
  using allowed_dest_tags = tmpl::list<>;
  template <typename Frame>
  using required_dest_tags = tmpl::list<>;
};

// Action that (artificially, for the test) clears the VolumeVarsInfo
template <typename TemporalIdTag>
struct ClearVolumeVarsInfo {
  template <
      typename ParallelComponent, typename DbTags, typename Metavariables,
      typename ArrayIndex,
      Requires<tmpl::list_contains_v<DbTags, intrp::Tags::NumberOfElements>> =
          nullptr>
  static void apply(db::DataBox<DbTags>& box,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/) {
    db::mutate<intrp::Tags::VolumeVarsInfo<Metavariables, TemporalIdTag>>(
        make_not_null(&box),
        [](const gsl::not_null<typename intrp::Tags::VolumeVarsInfo<
               Metavariables, TemporalIdTag>::type*>
               container) { container->clear(); });
  }
};

// Action that (artificially, for the test) inserts the given
// temporal_id into temporal_ids_when_data_has_been_interpolated.
template <typename InterpolationTargetTag>
struct AddToTemporalIdsWhenDataHasBeenInterpolated {
  template <
      typename ParallelComponent, typename DbTags, typename Metavariables,
      typename ArrayIndex,
      Requires<tmpl::list_contains_v<DbTags, intrp::Tags::NumberOfElements>> =
          nullptr>
  static void apply(
      db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/,
      const typename InterpolationTargetTag::temporal_id::type& temporal_id) {
    db::mutate<intrp::Tags::InterpolatedVarsHolders<Metavariables>>(
        make_not_null(&box),
        [&temporal_id](
            const gsl::not_null<typename intrp::Tags::InterpolatedVarsHolders<
                Metavariables>::type*>
                holders) {
          get<intrp::Vars::HolderTag<InterpolationTargetTag, Metavariables>>(
              *holders)
              .temporal_ids_when_data_has_been_interpolated.push_back(
                  temporal_id);
        });
  }
};

template <typename InterpolationTargetTag>
struct MockInterpolationTargetReceiveVars {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex, typename TemporalId,
            Requires<tmpl::list_contains_v<
                DbTags, intrp::Tags::TemporalIds<TemporalId>>> = nullptr>
  static void apply(
      db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/,
      const std::vector<::Variables<
          typename InterpolationTargetTag::vars_to_interpolate_to_target>>&
          vars_src,
      const std::vector<std::vector<size_t>>& global_offsets,
      const TemporalId& /*temporal_id*/) {
    size_t number_of_interpolated_points = 0;
    for (size_t i = 0; i < global_offsets.size(); ++i) {
      Scalar<DataVector> expected_vars{global_offsets[i].size()};
      number_of_interpolated_points += global_offsets[i].size();
      for (size_t s = 0; s < global_offsets[i].size(); ++s) {
        // Coords at this point. They are the same as the input coordinates,
        // but in strange order because of global_offsets.
        std::array<double, Metavariables::volume_dim> coords{
            {1.0 + 0.1 * global_offsets[i][s],
             1.0 + 0.12 * global_offsets[i][s],
             1.0 + 0.14 * global_offsets[i][s]}};
        const double lapse =
            2.0 * get<0>(coords) + 3.0 * get<1>(coords) +
            5.0 * get<2>(coords);  // Same formula as input lapse.
        get<>(expected_vars)[s] = square(lapse);
      }
      // We don't have that many points, so interpolation is good for
      // only a few digits.
      Approx custom_approx = Approx::custom().epsilon(1.e-5).scale(1.0);
      CHECK_ITERABLE_CUSTOM_APPROX(
          expected_vars, get<Tags::Square>(vars_src[i]), custom_approx);
    }
    // Make sure we have interpolated at the correct number of points.
    CHECK(number_of_interpolated_points == 15);
    // Change something in the DataBox so we can test that this function was
    // indeed called.  Put some unusual temporal_id into Tags::TemporalIds.
    // This is not the usual usage of Tags::TemporalIds; this is done just
    // for the test.
    Slab slab(0.0, 1.0);
    TimeStepId strange_temporal_id(true, 0, Time(slab, Rational(111, 135)));
    db::mutate<intrp::Tags::TemporalIds<TemporalId>>(
        make_not_null(&box),
        [&strange_temporal_id](
            const gsl::not_null<std::deque<TemporalId>*> temporal_ids) {
          temporal_ids->push_back(strange_temporal_id);
        });
  }
};

template <typename Metavariables, typename InterpolationTargetTag>
struct mock_interpolation_target {
  static_assert(
      tt::assert_conforms_to<InterpolationTargetTag,
                             intrp::protocols::InterpolationTargetTag>);
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = size_t;
  using component_being_mocked =
      intrp::InterpolationTarget<Metavariables, InterpolationTargetTag>;
  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Metavariables::volume_dim>>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<Actions::SetupDataBox,
                     intrp::Actions::InitializeInterpolationTarget<
                         Metavariables, InterpolationTargetTag>>>,
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::Registration, tmpl::list<>>,
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::Testing, tmpl::list<>>>;

  using replace_these_simple_actions =
      tmpl::list<intrp::Actions::InterpolationTargetReceiveVars<
          typename Metavariables::InterpolationTargetA>>;
  using with_these_simple_actions =
      tmpl::list<MockInterpolationTargetReceiveVars<
          typename Metavariables::InterpolationTargetA>>;
};

template <typename Metavariables>
struct mock_interpolator {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = size_t;
  using simple_tags = typename intrp::Actions::InitializeInterpolator<
      intrp::Tags::VolumeVarsInfo<Metavariables, ::Tags::TimeStepId>,
      intrp::Tags::InterpolatedVarsHolders<Metavariables>>::simple_tags;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>,
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::Registration, tmpl::list<>>,
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::Testing, tmpl::list<>>>;
  using component_being_mocked = void;  // not needed.
};

struct MockMetavariables {
  struct InterpolationTargetA
      : tt::ConformsTo<intrp::protocols::InterpolationTargetTag> {
    using temporal_id = ::Tags::TimeStepId;
    using vars_to_interpolate_to_target = tmpl::list<Tags::Square>;
    using compute_vars_to_interpolate = ComputeSquare;
    using compute_items_on_target = tmpl::list<>;
    using compute_target_points =
        ::intrp::TargetPoints::LineSegment<InterpolationTargetA, 3>;
    using post_interpolation_callback =
        intrp::callbacks::ObserveTimeSeriesOnSurface<tmpl::list<>,
                                                     InterpolationTargetA>;
  };
  using interpolator_source_vars = tmpl::list<gr::Tags::Lapse<DataVector>>;
  using interpolation_target_tags = tmpl::list<InterpolationTargetA>;
  static constexpr size_t volume_dim = 3;
  using component_list = tmpl::list<
      mock_interpolation_target<MockMetavariables, InterpolationTargetA>,
      mock_interpolator<MockMetavariables>>;
  enum class Phase { Initialization, Registration, Testing, Exit };
};

// Create volume data and send it to the interpolator.
template <typename interp_component, typename Metavariables,
          typename DomainCreatorType, typename DomainType>
void create_volume_data_and_send_it_to_interpolator(
    const gsl::not_null<ActionTesting::MockRuntimeSystem<Metavariables>*>
        runner,
    const DomainCreatorType& domain_creator, const DomainType& domain,
    const std::vector<ElementId<3>>& element_ids,
    const TimeStepId& temporal_id) {
  for (const auto& element_id : element_ids) {
    const auto& block = domain.blocks()[element_id.block_id()];
    ::Mesh<3> mesh{domain_creator.initial_extents()[element_id.block_id()],
                   Spectral::Basis::Legendre,
                   Spectral::Quadrature::GaussLobatto};
    if (block.is_time_dependent()) {
      ERROR("The block must be time-independent");
    }
    ElementMap<3, Frame::Inertial> map{element_id,
                                       block.stationary_map().get_clone()};
    const auto inertial_coords = map(logical_coordinates(mesh));
    ::Variables<typename Metavariables::interpolator_source_vars> output_vars(
        mesh.number_of_grid_points());
    auto& lapse = get<gr::Tags::Lapse<DataVector>>(output_vars);

    // Fill lapse with some analytic solution.
    get<>(lapse) = 2.0 * get<0>(inertial_coords) +
                   3.0 * get<1>(inertial_coords) +
                   5.0 * get<2>(inertial_coords);

    // Call the action on each element_id.
    runner->template simple_action<
        interp_component,
        ::intrp::Actions::InterpolatorReceiveVolumeData<
            typename Metavariables::InterpolationTargetA::temporal_id>>(
        0, temporal_id, element_id, mesh, std::move(output_vars));
  }
}

SPECTRE_TEST_CASE("Unit.NumericalAlgorithms.Interpolator.ReceiveVolumeData",
                  "[Unit]") {
  domain::creators::register_derived_with_charm();

  using metavars = MockMetavariables;
  using temporal_id_type =
      typename metavars::InterpolationTargetA::temporal_id::type;
  using target_component =
      mock_interpolation_target<metavars, metavars::InterpolationTargetA>;
  using interp_component = mock_interpolator<metavars>;

  // Make an InterpolatedVarsHolders containing the target points.
  const auto domain_creator =
      domain::creators::Shell(0.9, 4.9, 1, {{7, 7}}, false);
  const auto domain = domain_creator.create_domain();
  Slab slab(0.0, 1.0);
  TimeStepId temporal_id(true, 0, Time(slab, Rational(11, 15)));
  auto vars_holders = [&domain, &temporal_id]() {
    const size_t n_pts = 15;
    tnsr::I<DataVector, 3, Frame::Inertial> points(n_pts);
    for (size_t d = 0; d < 3; ++d) {
      for (size_t i = 0; i < n_pts; ++i) {
        points.get(d)[i] = 1.0 + (0.1 + 0.02 * d) * i;  // Chosen by hand.
      }
    }
    auto coords = block_logical_coordinates(domain, points);
    typename intrp::Tags::InterpolatedVarsHolders<metavars>::type
        vars_holders_l{};
    auto& vars_infos =
        get<intrp::Vars::HolderTag<metavars::InterpolationTargetA, metavars>>(
            vars_holders_l)
            .infos;
    vars_infos.emplace(std::make_pair(
        temporal_id,
        intrp::Vars::Info<3, typename metavars::InterpolationTargetA::
                                 vars_to_interpolate_to_target>{
            std::move(coords)}));
    return vars_holders_l;
  }();

  ActionTesting::MockRuntimeSystem<metavars> runner{
      {domain_creator.create_domain()}};
  ActionTesting::emplace_component_and_initialize<interp_component>(
      &runner, 0,
      {0_st,
       typename intrp::Tags::VolumeVarsInfo<
           metavars,
           typename metavars::InterpolationTargetA::temporal_id>::type{},
       typename intrp::Tags::InterpolatedVarsHolders<metavars>::type{
           vars_holders}});
  ActionTesting::emplace_component<target_component>(&runner, 0);
  for (size_t i = 0; i < 2; ++i) {
    ActionTesting::next_action<target_component>(make_not_null(&runner), 0);
  }
  ActionTesting::set_phase(make_not_null(&runner),
                           metavars::Phase::Registration);

  // Create Element_ids.
  std::vector<ElementId<3>> element_ids{};
  for (const auto& block : domain.blocks()) {
    const auto initial_ref_levs =
        domain_creator.initial_refinement_levels()[block.id()];
    auto elem_ids = initial_element_ids(block.id(), initial_ref_levs);
    element_ids.insert(element_ids.end(), elem_ids.begin(), elem_ids.end());
  }

  // Tell the interpolator how many elements there are by registering
  // each one.
  for (size_t i = 0; i < element_ids.size(); ++i) {
    runner.simple_action<interp_component, intrp::Actions::RegisterElement>(0);
  }
  ActionTesting::set_phase(make_not_null(&runner), metavars::Phase::Testing);

  create_volume_data_and_send_it_to_interpolator<interp_component>(
      make_not_null(&runner), domain_creator, domain, element_ids, temporal_id);

  // Should be no temporal_ids in the target box, since we never
  // put any there.
  CHECK(ActionTesting::get_databox_tag<
            target_component, intrp::Tags::TemporalIds<temporal_id_type>>(
            runner, 0)
            .empty());

  // Should be one queued simple action, MockInterpolationTargetReceiveVars.
  runner.invoke_queued_simple_action<target_component>(0);

  // Make sure that MockInterpolationTargetReceiveVars was called,
  // by looking for a funny temporal_id that it inserts for the specific
  // purpose of this test.
  CHECK(ActionTesting::get_databox_tag<
            target_component, intrp::Tags::TemporalIds<temporal_id_type>>(
            runner, 0)
            .front() ==
        TimeStepId(true, 0, Time(Slab(0.0, 1.0), Rational(111, 135))));

  // No more queued simple actions.
  CHECK(runner.is_simple_action_queue_empty<target_component>(0));

  // VolumeVarsInfo should be full now, since we didn't actually
  // do interpolation for this test.
  const auto& volume_vars_info = ActionTesting::get_databox_tag<
      interp_component,
      intrp::Tags::VolumeVarsInfo<
          metavars, typename metavars::InterpolationTargetA::temporal_id>>(
      runner, 0);
  CHECK(volume_vars_info.size() == 1);
  CHECK(volume_vars_info.at(temporal_id).size() == element_ids.size());

  // Now we will test that if temporal_ids_when_data_has_been_interpolated
  // is set for this temporal_id, subsequent calls of
  // InterpolatorReceiveVolumeData have no effect.
  // Do do this, we
  // 1) artificially clear volume_vars_info
  // 2) artificially set temporal_ids_when_data_has_been_interpolated
  // 3) call InterpolatorReceiveVolumeData
  // 4) Check that volume_vars_info is still empty, and that no simple_actions
  //    are queued.
  //
  runner
      .simple_action<interp_component,
                     ClearVolumeVarsInfo<
                         typename metavars::InterpolationTargetA::temporal_id>>(
          0);
  runner.simple_action<interp_component,
                       AddToTemporalIdsWhenDataHasBeenInterpolated<
                           metavars::InterpolationTargetA>>(0, temporal_id);

  create_volume_data_and_send_it_to_interpolator<interp_component>(
      make_not_null(&runner), domain_creator, domain, element_ids, temporal_id);

  // volume_vars_info should still be empty.
  CHECK(volume_vars_info.empty());

  // No more queued simple actions.
  CHECK(runner.is_simple_action_queue_empty<target_component>(0));
}
}  // namespace
