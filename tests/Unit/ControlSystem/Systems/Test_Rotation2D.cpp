// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

#include "ControlSystem/ApparentHorizons/Measurements.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Rotation.hpp"
#include "ControlSystem/Systems/Rotation.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/UpdateControlSystem.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation3D.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/PointwiseFunctions/PostNewtonian/BinaryTrajectories.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/UpdateMessageQueue.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

#include <cmath>
#include <ostream>
#include "Parallel/Printf.hpp"
#include "Utilities/EqualWithinRoundoff.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace control_system {
namespace {
struct MockMetavars;

constexpr size_t rotation_deriv_order = 2;
bool newtonian = false;

template <size_t DerivOrder>
struct MockRotation2DSystem : Rotation<DerivOrder> {
  struct MeasurementQueue : db::SimpleTag {
    using type = LinkedMessageQueue<
        double, tmpl::list<QueueTags::Center<ah::HorizonLabel::AhA>,
                           QueueTags::Center<ah::HorizonLabel::AhB>>>;
  };

  static constexpr size_t deriv_order_controlled = DerivOrder;

  struct process_measurement {
    template <typename Component, ah::HorizonLabel Horizon,
              typename Metavariables>
    static void apply(
        gsl::not_null<ActionTesting::MockRuntimeSystem<Metavariables>*> runner,
        const std::array<double, 3>& center_array,
        const LinkedMessageId<double>& measurement_id) {
      const DataVector center = array_to_datavector(center_array);

      ActionTesting::simple_action<
          Component, ::Actions::UpdateMessageQueue<
                         QueueTags::Center<Horizon>, MeasurementQueue,
                         UpdateControlSystem<deriv_order_controlled,
                                             ControlErrors::Rotation2D>>>(
          runner, 0, measurement_id, center);
    }
  };
};

template <typename ControlSystem>
using init_simple_tags = tmpl::list<
    control_system::Tags::Averager<ControlSystem::deriv_order_controlled - 1>,
    control_system::Tags::TimescaleTuner,
    control_system::Tags::ControlSystemName,
    control_system::Tags::Controller<ControlSystem::deriv_order_controlled>,
    typename ControlSystem::MeasurementQueue>;

template <typename Metavariables, typename MockControlSystem>
struct MockControlComponent {
  using array_index = int;
  using component_being_mocked =
      ControlComponent<Metavariables, MockControlSystem>;
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockSingletonChare;

  using simple_tags = init_simple_tags<MockControlSystem>;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      typename metavariables::Phase, metavariables::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>>;
};

template <typename Metavars>
struct MockElementComponent {
  using array_index = int;
  using chare_type = ActionTesting::MockArrayChare;

  using metavariables = Metavars;

  using simple_tags = tmpl::list<>;

  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Metavars::volume_dim>>;

  using mutable_global_cache_tags = tmpl::list<domain::Tags::FunctionsOfTime>;

  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename metavariables::Phase,
                                        metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};

template <typename Metavars>
using rotation_2d_component =
    MockControlComponent<Metavars, MockRotation2DSystem<rotation_deriv_order>>;
template <typename Metavars>
using element_component = MockElementComponent<Metavars>;

struct MockMetavars {
  static constexpr size_t volume_dim = 3;

  enum class Phase { Initialization, Testing, Exit };

  using component_list = tmpl::list<rotation_2d_component<MockMetavars>,
                                    element_component<MockMetavars>>;
};

using RotationMap = domain::CoordinateMaps::TimeDependent::Rotation<2>;
using Identity = domain::CoordinateMaps::Identity<1>;
using RotationMap3D =
    domain::CoordinateMaps::TimeDependent::ProductOf2Maps<RotationMap,
                                                          Identity>;

using CoordMap =
    domain::CoordinateMap<Frame::Grid, Frame::Inertial, RotationMap3D>;

std::array<double, 3> abs(const std::array<double, 3>& a) {
  return {std::abs(a[0]), std::abs(a[1]), std::abs(a[2])};
}

SPECTRE_TEST_CASE("Unit.ControlSystem.Systems.Rotation2D",
                  "[ControlSystem][Unit]") {
  // Global things
  domain::FunctionsOfTime::register_derived_with_charm();
  const double initial_separation = 15.0;
  double time = 0.0;
  std::optional<double> prev_time{};
  const double dt = 0.01;
  const double final_time = 500.0;

  // We don't need a real domain, just one that has the correct excision
  // spheres because the control errors use the `excision_spheres()` member
  // of a domain. The names are chosen to match the BinaryCompactObject
  // domain, which the control errors were based on and have these specific
  // names hard-coded into them.
  Domain<3> fake_domain{
      {},
      {},
      {{"ObjectAExcisionSphere",
        ExcisionSphere<3>{1.0, {{-0.5 * initial_separation, 0.0, 0.0}}}},
       {"ObjectBExcisionSphere",
        ExcisionSphere<3>{1.0, {{+0.5 * initial_separation, 0.0, 0.0}}}}}};

  // Analytic orbits
  const std::array<double, 3> initial_velocity{0.0, -0.0, 0.0};
  const BinaryTrajectories binary_trajectories{initial_separation,
                                               initial_velocity, newtonian};
  if (newtonian) {
    Parallel::printf("Newtonian," + get_output(rotation_deriv_order) + "\n");
  } else {
    Parallel::printf("PN," + get_output(rotation_deriv_order) + "\n");
  }

  // Names are hard coded here because we aren't using the full parallel
  // component infrastructure to create the control systems (too much overhead).
  const std::string rotation_name = "Rotation";

  // Inputs for control systems. We need different controllers, averagers, and
  // timescale tuners because different initial parameters get set for each
  // before they are put into the DataBoxes
  // deriv_order-1 because for controlling an nth derivative, we only need n-1
  // averaged derivatives.
  const double damping_timescale = 2.0;
  // args: timescales, max_timescale, min_timescale,
  //       decrease_threshold, increase_threshold,
  //       increase_factor, decrease_factor
  TimescaleTuner rotation_tuner{
      {damping_timescale}, 10.0, 0.1, 4.0e-2, 1.0e-2, 1.01, 0.98};
  // This +1 because the rotation control system is special. We are
  // controlling omega, but the control error that is calculated is actually an
  // infinitesimal angle that the binary has rotated through. Thus the nth
  // derivative of omega is the n+1st derivative of this angle
  const double update_fraction = 0.3;
  Controller<rotation_deriv_order> rotation_controller{update_fraction};
  // Same reason the rotation averager is 1 greater than the
  // translation/expansion averager
  Averager<rotation_deriv_order - 1> rotation_averager{0.25, false};
  LinkedMessageQueue<double,
                     tmpl::list<QueueTags::Center<ah::HorizonLabel::AhA>,
                                QueueTags::Center<ah::HorizonLabel::AhB>>>
      empty_queue{};

  // Since we aren't using the control system initialization actions, we need to
  // do all the initialization stuff manually.
  const double rotation_min_timescale = min(rotation_tuner.current_timescale());
  rotation_controller.assign_time_between_updates(rotation_min_timescale);
  rotation_controller.set_initial_time(time);

  // To create the DataBoxes for each control component. The LinkedMessageQueue
  // doesn't need initial parameters so just use default contruction
  tuples::tagged_tuple_from_typelist<
      init_simple_tags<MockRotation2DSystem<rotation_deriv_order>>>
      init_rot_tuple{rotation_averager, rotation_tuner, rotation_name,
                     rotation_controller, empty_queue};

  // Initial parameters needed. Expiration times would normally be set during
  // option parsing, so we have to do them manually here instead.
  const double initial_rotation_expiration_time =
      time + rotation_controller.get_update_fraction() * rotation_min_timescale;
  const double initial_omega = 0.01;
  const std::array<DataVector, 1> init_quaternion{{{1.0, 0.0, 0.0, 0.0}}};
  auto init_func_rotation =
      make_array<rotation_deriv_order + 1, DataVector>(DataVector{1, 0.0});
  init_func_rotation[1] = DataVector{{initial_omega}};

  // Create functions of time
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      initial_functions_of_time{};
  initial_functions_of_time[rotation_name] = std::make_unique<
      domain::FunctionsOfTime::PiecewisePolynomial<rotation_deriv_order>>(
      time, init_func_rotation, initial_rotation_expiration_time);

  // Setup runner and all components
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<MockMetavars>;
  MockRuntimeSystem runner{{std::move(fake_domain)},
                           {std::move(initial_functions_of_time)}};
  ActionTesting::emplace_singleton_component_and_initialize<
      rotation_2d_component<MockMetavars>>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, init_rot_tuple);
  ActionTesting::emplace_array_component<element_component<MockMetavars>>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0);

  ActionTesting::set_phase<MockMetavars>(make_not_null(&runner),
                                         MockMetavars::Phase::Testing);

  // Create coordinate maps for mapping the PN trajectories to the "grid" frame
  // where the control system does its calculations
  RotationMap rotation_map{rotation_name};
  Identity identity{};
  RotationMap3D rotation_map_3d{rotation_map, identity};

  CoordMap coord_map{rotation_map_3d};

  // Allocate now and just reuse as we go
  LinkedMessageId<double> measurement_id{};
  std::pair<std::array<double, 3>, std::array<double, 3>> positions{};
  tnsr::I<double, 3, Frame::Inertial> inertial_position_of_a{};
  tnsr::I<double, 3, Frame::Inertial> inertial_position_of_b{};
  tnsr::I<double, 3, Frame::Grid> grid_position_of_a_tnsr{};
  tnsr::I<double, 3, Frame::Grid> grid_position_of_b_tnsr{};
  std::array<double, 3> grid_position_of_a{};
  std::array<double, 3> grid_position_of_b{};

  // Get the functions of time from the cache to use in the maps
  const auto& cache =
      ActionTesting::cache<element_component<MockMetavars>>(runner, 0);
  const auto& functions_of_time =
      Parallel::get<domain::Tags::FunctionsOfTime>(cache);

  // Start loop
  while (time < final_time) {
    // Setup the measurement Id. This would normally be created in the control
    // system event.
    measurement_id = LinkedMessageId<double>{time, prev_time};

    // This whole switching between tensors and arrays is annoying and
    // clunky, but it's the best that could be done at the moment without
    // completely changing BinaryTrajectories, so we'll just go with it for now.

    // Get trajectory in "inertial coordinates" as arrays
    positions = binary_trajectories.positions_no_expansion(time);

    // Covert arrays to tensor so we can pass them into the coordinate map
    for (size_t i = 0; i < 3; i++) {
      inertial_position_of_a.get(i) = gsl::at(positions.first, i);
      inertial_position_of_b.get(i) = gsl::at(positions.second, i);
    }

    // Convert to "grid coordinaes"
    grid_position_of_a_tnsr =
        *coord_map.inverse(inertial_position_of_a, time, functions_of_time);
    grid_position_of_b_tnsr =
        *coord_map.inverse(inertial_position_of_b, time, functions_of_time);

    // Convert tensors back to arrays so we can pass them to the control systems
    for (size_t i = 0; i < 3; i++) {
      gsl::at(grid_position_of_a, i) = grid_position_of_a_tnsr.get(i);
      gsl::at(grid_position_of_b, i) = grid_position_of_b_tnsr.get(i);
    }

    // Apply measurements for A
    MockRotation2DSystem<rotation_deriv_order>::process_measurement::apply<
        rotation_2d_component<MockMetavars>, ah::HorizonLabel::AhA>(
        make_not_null(&runner), grid_position_of_a, measurement_id);
    // Apply measurements for B
    MockRotation2DSystem<rotation_deriv_order>::process_measurement::apply<
        rotation_2d_component<MockMetavars>, ah::HorizonLabel::AhB>(
        make_not_null(&runner), grid_position_of_b, measurement_id);

    // At this point, the control systems for each transformation should have
    // done their thing and updated the functions of time (if they had enough
    // data).

    // Go to next time
    prev_time = time;
    time += dt;
  }

  // Get analytic position in inertial coordinates
  positions = binary_trajectories.positions_no_expansion(time);

  // Get position of objects in grid coordinates using the coordinate map that
  // has had its functions of time updated by the control system
  for (size_t i = 0; i < 3; i++) {
    inertial_position_of_a.get(i) = gsl::at(positions.first, i);
    inertial_position_of_b.get(i) = gsl::at(positions.second, i);
  }
  grid_position_of_a_tnsr =
      *coord_map.inverse(inertial_position_of_a, time, functions_of_time);
  grid_position_of_b_tnsr =
      *coord_map.inverse(inertial_position_of_b, time, functions_of_time);
  for (size_t i = 0; i < 3; i++) {
    gsl::at(grid_position_of_a, i) = grid_position_of_a_tnsr.get(i);
    gsl::at(grid_position_of_b, i) = grid_position_of_b_tnsr.get(i);
  }
  const std::array<double, 3> expected_grid_position_of_a{
      {-0.5 * initial_separation, 0.0, 0.0}};
  const std::array<double, 3> expected_grid_position_of_b{
      {0.5 * initial_separation, 0.0, 0.0}};

  const auto& updated_rotation_tuner =
      ActionTesting::get_databox_tag<rotation_2d_component<MockMetavars>,
                                     control_system::Tags::TimescaleTuner>(
          runner, 0);

  const auto& rotation_f_of_t = dynamic_cast<
      domain::FunctionsOfTime::PiecewisePolynomial<rotation_deriv_order>&>(
      *functions_of_time.at(rotation_name));

  std::ostringstream os;
  os << "Abs(expected - actual) at time " << final_time << ": "
     << abs(expected_grid_position_of_b - grid_position_of_b) << "\n";
  std::string f_of_t_str = get_output(rotation_f_of_t);
  os << "Number of times updated: "
     << std::count(f_of_t_str.begin(), f_of_t_str.end(), '=') / 2.0 << "\n";
  os << "Rotation timescales: " << updated_rotation_tuner.current_timescale()
     << "\n";
  os << get_output(rotation_f_of_t) << "\n";
  // os << get_output(expansion_f_of_t) << "\n";
  Parallel::printf(os.str());

  const auto omega = rotation_f_of_t.func_and_deriv(time)[1];
  const auto pn_omega = binary_trajectories.omega(time);
  // The control system gets more accurate the longer you run for.
  // This is the accuracy we can achieve in this amount of time.
  Approx custom_approx = Approx::custom().epsilon(5.0e-5).scale(1.0);
  const DataVector expected_omega{{pn_omega}};
  CHECK_ITERABLE_CUSTOM_APPROX(expected_omega, omega, custom_approx);

  CHECK_ITERABLE_CUSTOM_APPROX(expected_grid_position_of_a, grid_position_of_a,
                               custom_approx);
  CHECK_ITERABLE_CUSTOM_APPROX(expected_grid_position_of_b, grid_position_of_b,
                               custom_approx);
}
}  // namespace
}  // namespace control_system
