// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>

#include "ControlSystem/ApparentHorizons/Measurements.hpp"
#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Component.hpp"
#include "ControlSystem/ControlErrors/Expansion.hpp"
#include "ControlSystem/ControlErrors/Rotation.hpp"
#include "ControlSystem/ControlErrors/Translation.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/DataVectorHelpers.hpp"
#include "ControlSystem/Systems/Expansion.hpp"
#include "ControlSystem/Systems/Rotation.hpp"
#include "ControlSystem/Systems/Translation.hpp"
#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "ControlSystem/UpdateControlSystem.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/LinkedMessageQueue.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Domain.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Domain/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/TestingFramework.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system::TestHelpers {
template <typename ControlSystem>
using init_simple_tags =
    tmpl::list<control_system::Tags::Averager<ControlSystem::deriv_order - 1>,
               control_system::Tags::TimescaleTuner,
               control_system::Tags::ControlSystemName,
               control_system::Tags::Controller<ControlSystem::deriv_order>,
               control_system::Tags::ControlError<ControlSystem>,
               typename ControlSystem::MeasurementQueue>;

template <typename Metavariables, typename ControlSystem>
struct MockControlComponent {
  using array_index = int;
  using component_being_mocked = ControlComponent<Metavariables, ControlSystem>;
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockSingletonChare;

  using system = ControlSystem;

  using simple_tags = init_simple_tags<ControlSystem>;

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

  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize,
                 control_system::Tags::MeasurementTimescales>;

  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename metavariables::Phase,
                                        metavariables::Phase::Initialization,
                                        tmpl::list<>>>;
};

template <size_t TranslationDerivOrder, size_t RotationDerivOrder,
          size_t ExpansionDerivOrder>
struct MockMetavars {
  static constexpr size_t volume_dim = 3;

  enum class Phase { Initialization, Testing, Exit };

  using metavars = MockMetavars<TranslationDerivOrder, RotationDerivOrder,
                                ExpansionDerivOrder>;

  static constexpr bool using_translation = TranslationDerivOrder != 0;
  static constexpr bool using_rotation = RotationDerivOrder != 0;
  static constexpr bool using_expansion = ExpansionDerivOrder != 0;

  // Even if we aren't using certain control systems, we still need valid deriv
  // orders becuse everything is constructed by default in the SystemHelper. The
  // bools above just determine if the functions of time are actually created or
  // not because that's what matters
  static constexpr size_t trans_deriv_order =
      using_translation ? TranslationDerivOrder : 2;
  static constexpr size_t rot_deriv_order =
      using_rotation ? RotationDerivOrder : 2;
  static constexpr size_t exp_deriv_order =
      using_expansion ? ExpansionDerivOrder : 2;

  using translation_system =
      control_system::Systems::Translation<trans_deriv_order>;
  using rotation_system = control_system::Systems::Rotation<rot_deriv_order>;
  using expansion_system = control_system::Systems::Expansion<exp_deriv_order>;

  using element_component = MockElementComponent<metavars>;
  using translation_component =
      MockControlComponent<metavars, translation_system>;
  using rotation_component = MockControlComponent<metavars, rotation_system>;
  using expansion_component = MockControlComponent<metavars, expansion_system>;

  using control_components = tmpl::flatten<tmpl::list<
      tmpl::conditional_t<using_translation, tmpl::list<translation_component>,
                          tmpl::list<>>,
      tmpl::conditional_t<using_rotation, tmpl::list<rotation_component>,
                          tmpl::list<>>,
      tmpl::conditional_t<using_expansion, tmpl::list<expansion_component>,
                          tmpl::list<>>>>;

  using component_list =
      tmpl::flatten<tmpl::list<element_component, control_components>>;
};

/*!
 * \brief Helper struct for testing control systems
 *
 * To signify which control systems you want, set the corresponding
 * DerivOrder. To turn control systems off, put 0 for their DerivOrder in the
 * templates of the metavariables. For example, if you only want rotation
 * control with deriv order 2, you'd set up the metavariables as:
 * MockMetavars<0, 2, 0>{};
 *
 * Ideally we'd construct the runner here and just pass that to the test to
 * simplify as must of the work as possible, but MockRuntimeSystems aren't
 * copy- or move-able so we have to make the necessary info available. The
 * simplist way to do this was to have functions that return references to the
 * memeber variables.
 */
template <typename Metavars>
struct SystemHelper {
  static constexpr size_t trans_deriv_order = Metavars::trans_deriv_order;
  static constexpr size_t rot_deriv_order = Metavars::rot_deriv_order;
  static constexpr size_t exp_deriv_order = Metavars::exp_deriv_order;

  static constexpr bool using_translation = Metavars::using_translation;
  static constexpr bool using_rotation = Metavars::using_rotation;
  static constexpr bool using_expansion = Metavars::using_expansion;

  using translation_system = typename Metavars::translation_system;
  using rotation_system = typename Metavars::rotation_system;
  using expansion_system = typename Metavars::expansion_system;

  using element_component = typename Metavars::element_component;
  using control_components = typename Metavars::control_components;

  using translation_init_simple_tags = init_simple_tags<translation_system>;
  using rotation_init_simple_tags = init_simple_tags<rotation_system>;
  using expansion_init_simple_tags = init_simple_tags<expansion_system>;

  // Memebers that may be moved out of this struct once they are
  // constructed
  auto& domain() { return domain_; }
  auto& initial_functions_of_time() { return initial_functions_of_time_; }
  auto& initial_measurement_timescales() {
    return initial_measurement_timescales_;
  }

  // Memebers that won't be moved out of this struct
  const auto& init_trans_tuple() { return init_trans_tuple_; }
  const auto& init_rot_tuple() { return init_rot_tuple_; }
  const auto& init_exp_tuple() { return init_exp_tuple_; }
  const auto& grid_position_of_a() { return grid_position_of_a_; }
  const auto& grid_position_of_b() { return grid_position_of_b_; }
  const auto& translation_name() { return translation_name_; }
  const auto& rotation_name() { return rotation_name_; }
  const auto& expansion_name() { return expansion_name_; }

  void setup_control_system_test(const double initial_time,
                                 const double initial_separation) {
    domain::FunctionsOfTime::register_derived_with_charm();
    initial_time_ = initial_time;

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
    domain_ = std::move(fake_domain);

    // Inputs for control systems. We need different controllers, averagers, and
    // timescale tuners because different initial parameters get set for each
    // before they are put into the DataBoxes
    Averager<trans_deriv_order - 1> translation_averager{0.25, false};
    Averager<exp_deriv_order - 1> expansion_averager{0.25, false};
    Averager<rot_deriv_order - 1> rotation_averager{0.25, false};
    const double update_timescale = 0.5;
    TimescaleTuner translation_tuner{
        {update_timescale, update_timescale, update_timescale},
        10.0,
        0.1,
        4.0e-3,
        1.0e-3,
        1.01,
        0.98};
    TimescaleTuner rotation_tuner{
        {update_timescale, update_timescale, update_timescale},
        10.0,
        0.1,
        4.0e-2,
        1.0e-2,
        1.01,
        0.98};
    TimescaleTuner expansion_tuner{
        {update_timescale}, 10.0, 0.1, 2.0, 0.1, 1.01, 0.99};
    Controller<trans_deriv_order> translation_controller{0.3};
    Controller<rot_deriv_order> rotation_controller{0.3};
    Controller<exp_deriv_order> expansion_controller{0.3};

    ControlErrors::Translation translation_control_error{};
    ControlErrors::Rotation rotation_control_error{};
    ControlErrors::Expansion expansion_control_error{};

    // Since we aren't using the control system initialization actions, we need
    // to do all the initialization stuff manually.
    const double translation_min_timescale =
        min(translation_tuner.current_timescale());
    const double rotation_min_timescale =
        min(rotation_tuner.current_timescale());
    const double expansion_min_timescale =
        min(expansion_tuner.current_timescale());
    translation_controller.assign_time_between_updates(
        translation_min_timescale);
    translation_controller.set_initial_time(initial_time_);
    rotation_controller.assign_time_between_updates(rotation_min_timescale);
    rotation_controller.set_initial_time(initial_time_);
    expansion_controller.assign_time_between_updates(expansion_min_timescale);
    expansion_controller.set_initial_time(initial_time_);
    const std::array<DataVector, 1> translation_measurement_timescale{
        {control_system::Tags::detail::calculate_measurement_timescales(
            translation_controller, translation_tuner)}};
    const std::array<DataVector, 1> rotation_measurement_timescale{
        {control_system::Tags::detail::calculate_measurement_timescales(
            rotation_controller, rotation_tuner)}};
    const std::array<DataVector, 1> expansion_measurement_timescale{
        {control_system::Tags::detail::calculate_measurement_timescales(
            expansion_controller, expansion_tuner)}};
    translation_averager.assign_time_between_measurements(
        min(translation_measurement_timescale[0]));
    rotation_averager.assign_time_between_measurements(
        min(rotation_measurement_timescale[0]));
    expansion_averager.assign_time_between_measurements(
        min(expansion_measurement_timescale[0]));

    // To create the DataBoxes for each control component. The
    // LinkedMessageQueue doesn't need initial parameters so just use default
    // contruction
    init_trans_tuple_ =
        tuples::tagged_tuple_from_typelist<translation_init_simple_tags>{
            translation_averager,      translation_tuner,
            translation_name_,         translation_controller,
            translation_control_error, empty_queue_};
    init_rot_tuple_ =
        tuples::tagged_tuple_from_typelist<rotation_init_simple_tags>{
            rotation_averager,   rotation_tuner,         rotation_name_,
            rotation_controller, rotation_control_error, empty_queue_};
    init_exp_tuple_ =
        tuples::tagged_tuple_from_typelist<expansion_init_simple_tags>{
            expansion_averager,   expansion_tuner,         expansion_name_,
            expansion_controller, expansion_control_error, empty_queue_};

    // Initial parameters needed. Expiration times would normally be set during
    // option parsing, so we have to do them manually here instead.
    [[maybe_unused]] const double initial_translation_expiration_time =
        translation_controller.get_update_fraction() *
        translation_min_timescale;
    [[maybe_unused]] const double initial_rotation_expiration_time =
        rotation_controller.get_update_fraction() * rotation_min_timescale;
    [[maybe_unused]] const double initial_expansion_expiration_time =
        expansion_controller.get_update_fraction() * expansion_min_timescale;
    const double initial_expansion = 1.0;
    const double initial_omega = 0.01;
    const double expansion_velocity_outer_boundary = 0.0;
    const double decay_timescale_outer_boundary = 0.05;
    auto init_func_expansion =
        make_array<exp_deriv_order + 1, DataVector>(DataVector{1, 0.0});
    init_func_expansion[0][0] = initial_expansion;
    auto init_func_translation =
        make_array<trans_deriv_order + 1, DataVector>(DataVector{3, 0.0});
    const std::array<DataVector, 1> init_quaternion{{{1.0, 0.0, 0.0, 0.0}}};
    auto init_func_rotation =
        make_array<rot_deriv_order + 1, DataVector>(DataVector{3, 0.0});
    init_func_rotation[1] = DataVector{{0.0, 0.0, initial_omega}};

    if constexpr (using_translation) {
      initial_functions_of_time_[translation_name_] = std::make_unique<
          domain::FunctionsOfTime::PiecewisePolynomial<trans_deriv_order>>(
          initial_time_, init_func_translation,
          initial_translation_expiration_time);
      initial_measurement_timescales_[translation_name_] =
          std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
              initial_time_, translation_measurement_timescale,
              initial_translation_expiration_time);

      // Translation control error requires a quaternion and expansion thus we
      // have to add these in (otherwise we'd be adding in the whole control
      // systems themselves which we don't want). To avoid solving an ODE, we
      // use a constant PiecewisePolynomial representing the unit quaternion to
      // signify there isn't any rotation
      auto local_init_func_rotation =
          make_array<1, DataVector>(DataVector{4, 0.0});
      auto local_init_func_expansion =
          make_array<1, DataVector>(DataVector{1, 0.0});
      local_init_func_rotation[0][0] = 1.0;
      local_init_func_expansion[0][0] = 1.0;
      initial_functions_of_time_[rotation_name_] =
          std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
              initial_time, local_init_func_rotation,
              std::numeric_limits<double>::infinity());
      initial_functions_of_time_[expansion_name_] =
          std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
              initial_time, local_init_func_expansion,
              std::numeric_limits<double>::infinity());
    }
    if constexpr (using_rotation) {
      initial_functions_of_time_[rotation_name_] = std::make_unique<
          domain::FunctionsOfTime::QuaternionFunctionOfTime<rot_deriv_order>>(
          initial_time_, init_quaternion, init_func_rotation,
          initial_rotation_expiration_time);
      initial_measurement_timescales_[rotation_name_] =
          std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
              initial_time_, rotation_measurement_timescale,
              initial_rotation_expiration_time);
    }
    if constexpr (using_expansion) {
      initial_functions_of_time_[expansion_name_] = std::make_unique<
          domain::FunctionsOfTime::PiecewisePolynomial<exp_deriv_order>>(
          initial_time_, init_func_expansion,
          initial_expansion_expiration_time);
      initial_functions_of_time_[expansion_name_ + "OuterBoundary"s] =
          std::make_unique<domain::FunctionsOfTime::FixedSpeedCubic>(
              initial_expansion, initial_time_,
              expansion_velocity_outer_boundary,
              decay_timescale_outer_boundary);
      initial_measurement_timescales_[expansion_name_] =
          std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<0>>(
              initial_time_, expansion_measurement_timescale,
              initial_expansion_expiration_time);
    }
  }

  template <typename Generator, typename F, typename CoordMap>
  void run_control_system_test(
      ActionTesting::MockRuntimeSystem<Metavars>& runner,
      const double final_time, gsl::not_null<Generator*> generator,
      const F position_function, const CoordMap& coord_map) {
    // Allocate now and just reuse as we go
    LinkedMessageId<double> measurement_id{};
    Strahlkorper<Frame::Grid> strahlkorper_a{};
    Strahlkorper<Frame::Grid> strahlkorper_b{};
    std::pair<std::array<double, 3>, std::array<double, 3>> positions{};
    tnsr::I<double, 3, Frame::Inertial> inertial_position_of_a{};
    tnsr::I<double, 3, Frame::Inertial> inertial_position_of_b{};
    tnsr::I<double, 3, Frame::Grid> grid_position_of_a_tnsr{};
    tnsr::I<double, 3, Frame::Grid> grid_position_of_b_tnsr{};
    double time = initial_time_;
    std::optional<double> prev_time{};
    double dt = 0.0;

    auto& cache = ActionTesting::cache<element_component>(runner, 0);
    const auto& functions_of_time =
        Parallel::get<domain::Tags::FunctionsOfTime>(cache);
    const auto& measurement_timescales =
        Parallel::get<control_system::Tags::MeasurementTimescales>(cache);

    // Start loop
    while (time < final_time) {
      // Setup the measurement Id. This would normally be created in the control
      // system event.
      measurement_id = LinkedMessageId<double>{time, prev_time};

      // This whole switching between tensors and arrays is annoying and
      // clunky, but it's the best that could be done at the moment without
      // completely changing BinaryTrajectories, so we'll just go with it for
      // now.

      // Get trajectory in "inertial coordinates" as arrays
      positions = position_function(time);

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

      // Convert tensors back to arrays so we can pass them to the control
      // systems
      for (size_t i = 0; i < 3; i++) {
        gsl::at(grid_position_of_a_, i) = grid_position_of_a_tnsr.get(i);
        gsl::at(grid_position_of_b_, i) = grid_position_of_b_tnsr.get(i);
      }

      // Construct strahlkorpers to pass to control systems. Only the centers
      // matter.
      strahlkorper_a =
          Strahlkorper<Frame::Grid>{2, 2, 1.0, grid_position_of_a_};
      strahlkorper_b =
          Strahlkorper<Frame::Grid>{2, 2, 1.0, grid_position_of_b_};

      // Apply measurements
      tmpl::for_each<control_components>([&runner, &generator, &measurement_id,
                                          &strahlkorper_a, &cache,
                                          &strahlkorper_b](
                                             auto control_component) {
        using component = tmpl::type_from<decltype(control_component)>;
        using system = typename component::system;
        system::process_measurement::apply(
            ah::BothHorizons::FindHorizon<ah::HorizonLabel::AhA>{},
            strahlkorper_a, cache, measurement_id);
        CHECK(ActionTesting::number_of_queued_simple_actions<component>(
                  runner, 0) == 1);
        system::process_measurement::apply(
            ah::BothHorizons::FindHorizon<ah::HorizonLabel::AhB>{},
            strahlkorper_b, cache, measurement_id);
        CHECK(ActionTesting::number_of_queued_simple_actions<component>(
                  runner, 0) == 2);
        // We invoke a random measurement becuase during a normal simulation
        // we don't know which measurement will reach the control system
        // first because of charm++ communication
        ActionTesting::invoke_random_queued_simple_action<control_components>(
            make_not_null(&runner), generator,
            ActionTesting::array_indices_with_queued_simple_actions<
                control_components>(make_not_null(&runner)));
        ActionTesting::invoke_queued_simple_action<component>(
            make_not_null(&runner), 0);
      });

      // At this point, the control systems for each transformation should have
      // done their thing and updated the functions of time (if they had enough
      // data).

      // Our dt is set by the smallest measurement timescale. The control system
      // updates these timescales when it updates the functions of time
      prev_time = time;
      dt = std::numeric_limits<double>::max();
      for (auto& [name, measurement_timescale] : measurement_timescales) {
        dt = std::min(dt, min(measurement_timescale->func(time)[0]));
      }
      time += dt;
    }

    // Get analytic position in inertial coordinates
    positions = position_function(final_time);

    // Get position of objects in grid coordinates using the coordinate map that
    // has had its functions of time updated by the control system
    for (size_t i = 0; i < 3; i++) {
      inertial_position_of_a.get(i) = gsl::at(positions.first, i);
      inertial_position_of_b.get(i) = gsl::at(positions.second, i);
    }
    grid_position_of_a_tnsr = *coord_map.inverse(inertial_position_of_a,
                                                 final_time, functions_of_time);
    grid_position_of_b_tnsr = *coord_map.inverse(inertial_position_of_b,
                                                 final_time, functions_of_time);
    for (size_t i = 0; i < 3; i++) {
      gsl::at(grid_position_of_a_, i) = grid_position_of_a_tnsr.get(i);
      gsl::at(grid_position_of_b_, i) = grid_position_of_b_tnsr.get(i);
    }
  }

 private:
  // Memebers that may be moved out of this struct once they are
  // constructed
  Domain<3> domain_;
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      initial_functions_of_time_{};
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      initial_measurement_timescales_{};

  // Memebers that won't be moved out of this struct
  LinkedMessageQueue<double,
                     tmpl::list<QueueTags::Center<ah::HorizonLabel::AhA>,
                                QueueTags::Center<ah::HorizonLabel::AhB>>>
      empty_queue_{};
  tuples::tagged_tuple_from_typelist<translation_init_simple_tags>
      init_trans_tuple_{Averager<trans_deriv_order - 1>{},
                        TimescaleTuner{},
                        "",
                        Controller<trans_deriv_order>{},
                        ControlErrors::Translation{},
                        empty_queue_};
  tuples::tagged_tuple_from_typelist<rotation_init_simple_tags> init_rot_tuple_{
      Averager<rot_deriv_order - 1>{}, TimescaleTuner{},          "",
      Controller<rot_deriv_order>{},   ControlErrors::Rotation{}, empty_queue_};
  tuples::tagged_tuple_from_typelist<expansion_init_simple_tags>
      init_exp_tuple_{Averager<exp_deriv_order - 1>{},
                      TimescaleTuner{},
                      "",
                      Controller<exp_deriv_order>{},
                      ControlErrors::Expansion{},
                      empty_queue_};
  std::array<double, 3> grid_position_of_a_{};
  std::array<double, 3> grid_position_of_b_{};
  const std::string translation_name_{translation_system::name()};
  const std::string rotation_name_{rotation_system::name()};
  const std::string expansion_name_{expansion_system::name()};
  double initial_time_{std::numeric_limits<double>::signaling_NaN()};
};
}  // namespace control_system::TestHelpers
