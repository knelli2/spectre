// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "ControlSystem/Tags.hpp"
#include "ControlSystem/Tags/MeasurementTimescales.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/Identity.hpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Translation.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/ControlSystem/SystemHelpers.hpp"
#include "Helpers/PointwiseFunctions/PostNewtonian/BinaryTrajectories.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "NumericalAlgorithms/Interpolation/CubicSpline.hpp"
#include "Parallel/Phase.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StdArrayHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace control_system {
namespace {

using IdentityMap = domain::CoordinateMaps::Identity<1>;
using RotationMap2D = domain::CoordinateMaps::TimeDependent::Rotation<2>;
using RotationMap =
    domain::CoordinateMaps::TimeDependent::ProductOf2Maps<RotationMap2D,
                                                          IdentityMap>;
using ExpansionMap = domain::CoordinateMaps::TimeDependent::CubicScale<3>;

using CoordMap = domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                       ExpansionMap, RotationMap>;

std::string create_input_string(const std::string& name) {
  const std::string base_string1{"  "s + name +
                                 ":\n"
                                 "    Averager:\n"
                                 "      AverageTimescaleFraction: 0.25\n"
                                 "      Average0thDeriv: false\n"
                                 "    Controller:\n"
                                 "      UpdateFraction: 0.3\n"
                                 "    TimescaleTuner:\n"};
  const std::string base_string2{
      "      MinTimescale: 0.01\n"
      "      MaxTimescale: 20.0\n"
      "      DecreaseThreshold: 0.001\n"
      "      IncreaseThreshold: 0.00025\n"
      "      IncreaseFactor: 1.01\n"
      "      DecreaseFactor: 0.98\n"
      "    ControlError:\n"};

  const std::string rot_part =
      name == "Rotation"s ? "      OnlyAboutZAxis: true\n" : "";

  const std::string timescales =
      (name == "Translation"s) ? "      InitialTimescales: [0.2, 0.2, 0.2]\n"
                               : "      InitialTimescales: [0.2]\n";

  return base_string1 + timescales + base_string2 + rot_part;
}

std::pair<std::array<intrp::CubicSpline, 3>, std::array<intrp::CubicSpline, 3>>
read_spec_horizons(const std::string& filename, const std::string& /*prefix*/) {
  h5::H5File<h5::AccessType::ReadOnly> h5file{filename};
  const auto& aha_dat =
      h5file.get<h5::Dat>("/ApparentHorizons/ControlSystemAhA_Centers_sorted");

  const size_t num_rows = aha_dat.get_dimensions()[0];

  std::array<std::vector<double>, 3> aha_data{std::vector<double>(num_rows),
                                              std::vector<double>(num_rows),
                                              std::vector<double>(num_rows)};
  std::array<std::vector<double>, 3> ahb_data{std::vector<double>(num_rows),
                                              std::vector<double>(num_rows),
                                              std::vector<double>(num_rows)};
  std::vector<double> times(num_rows);

  const Matrix aha_all_data = aha_dat.get_data();
  h5file.close_current_object();

  const auto& ahb_dat =
      h5file.get<h5::Dat>("/ApparentHorizons/ControlSystemAhB_Centers_sorted");
  const Matrix ahb_all_data = ahb_dat.get_data();
  h5file.close_current_object();

  const size_t inertial_x_idx = 4;
  for (size_t i = 0; i < num_rows; i++) {
    times[i] = aha_all_data(i, 0);

    for (size_t j = 0; j < 3; j++) {
      gsl::at(aha_data, j)[i] = aha_all_data(i, inertial_x_idx + j);
      gsl::at(ahb_data, j)[i] = ahb_all_data(i, inertial_x_idx + j);
    }
  }

  std::array<intrp::CubicSpline, 3> aha_spline{};
  std::array<intrp::CubicSpline, 3> ahb_spline{};

  for (size_t i = 0; i < 3; i++) {
    aha_spline[i] = intrp::CubicSpline(times, aha_data[i]);
    ahb_spline[i] = intrp::CubicSpline(times, ahb_data[i]);
  }

  // Swap aha and ahb because SpEC does it backwards
  // But these are spectre so don't swap
  return std::make_pair(std::move(aha_spline), std::move(ahb_spline));
}

template <size_t TranslationDerivOrder, size_t RotationDerivOrder,
          size_t ExpansionDerivOrder>
void test_rotscaletrans_control_system() {
  INFO("Translation: "s + get_output(TranslationDerivOrder) + ", Rotation: "s +
       get_output(RotationDerivOrder) + ", Expansion: "s +
       get_output(ExpansionDerivOrder));
  using metavars =
      TestHelpers::MockMetavars<TranslationDerivOrder, RotationDerivOrder,
                                ExpansionDerivOrder>;
  using element_component = typename metavars::element_component;
  using rotation_component = typename metavars::rotation_component;
  using expansion_component = typename metavars::expansion_component;
  using rotation_system = typename metavars::rotation_system;
  using expansion_system = typename metavars::expansion_system;
  using obs_writer = typename metavars::observer_component;
  MAKE_GENERATOR(gen);

  // Global things
  domain::FunctionsOfTime::register_derived_with_charm();
  const double initial_time = 0.0;
  const double initial_separation = 15.366;
  // This final time is chosen so that the damping timescales have adequate time
  // to reach the maximum damping timescale
  const double final_time = 1.0;

  // Set up the system helper
  control_system::TestHelpers::SystemHelper<metavars> system_helper{};

  const std::string rotation_name =
      system_helper.template name<rotation_system>();
  const std::string expansion_name =
      system_helper.template name<expansion_system>();

  std::string input_options =
      "Evolution:\n"
      "  InitialTime: 0.0\n"
      "ControlSystems:\n"
      "  WriteDataToDisk: true\n"
      "  MeasurementsPerUpdate: 4\n";
  input_options += create_input_string(rotation_name);
  input_options += create_input_string(expansion_name);

  const auto initialize_functions_of_time =
      [](const gsl::not_null<std::unordered_map<
             std::string,
             std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>*>
             functions_of_time,
         const double local_initial_time,
         const std::unordered_map<std::string, double>&
             initial_expiration_times) {
        TestHelpers::initialize_expansion_functions_of_time<expansion_system>(
            functions_of_time, local_initial_time, initial_expiration_times);
        return TestHelpers::initialize_rotation2d_functions_of_time<
            rotation_system>(functions_of_time, local_initial_time,
                             initial_expiration_times);
      };

  // Initialize everything within the system helper
  system_helper.setup_control_system_test(initial_time, initial_separation,
                                          input_options,
                                          initialize_functions_of_time);

  // Get references to everything that was set up inside the system helper. The
  // domain and two functions of time are not const references because they need
  // to be moved into the runner
  auto& domain = system_helper.domain();
  auto& initial_functions_of_time = system_helper.initial_functions_of_time();
  auto& initial_measurement_timescales =
      system_helper.initial_measurement_timescales();
  const auto& init_rot_tuple =
      system_helper.template init_tuple<rotation_system>();
  const auto& init_exp_tuple =
      system_helper.template init_tuple<expansion_system>();

  // Setup runner and all components
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavars>;
  MockRuntimeSystem runner{{"DummyFileName", std::move(domain), 4},
                           {std::move(initial_functions_of_time),
                            std::move(initial_measurement_timescales)}};
  ActionTesting::emplace_singleton_component_and_initialize<rotation_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, init_rot_tuple);
  ActionTesting::emplace_singleton_component_and_initialize<
      expansion_component>(make_not_null(&runner), ActionTesting::NodeId{0},
                           ActionTesting::LocalCoreId{0}, init_exp_tuple);
  ActionTesting::emplace_array_component<element_component>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0);
  ActionTesting::emplace_nodegroup_component_and_initialize<obs_writer>(
      make_not_null(&runner), {});

  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  const std::string spec_filename{
      "/home/knelli/Documents/research/sims/mirror_spec/"
      "sorted_run_spectre_out.h5"};

  const auto horizons = read_spec_horizons(spec_filename, "Trajectory_"s);

  IdentityMap identity{};
  RotationMap2D rotation_2d{rotation_name};
  RotationMap rotation_map{rotation_2d, identity};
  ExpansionMap expansion_map{300.0, expansion_name,
                             expansion_name + "OuterBoundary"s};
  CoordMap coord_map{expansion_map, rotation_map};

  const auto position_function = [&horizons](const double time) {
    const auto& aha_spline = horizons.first;
    const auto& ahb_spline = horizons.second;

    std::array<double, 3> aha_position{aha_spline[0](time), aha_spline[1](time),
                                       aha_spline[2](time)};
    std::array<double, 3> ahb_position{ahb_spline[0](time), ahb_spline[1](time),
                                       ahb_spline[2](time)};

    return std::make_pair(std::move(aha_position), std::move(ahb_position));
  };

  const auto horizon_function = [&position_function, &coord_map,
                                 &runner](const double time) {
    // const auto& aha_spline = horizons.first;
    // const auto& ahb_spline = horizons.second;

    // std::array<double, 3> aha_position{aha_spline[0](time),
    // aha_spline[1](time),
    //                                    aha_spline[2](time)};
    // std::array<double, 3> ahb_position{ahb_spline[0](time),
    // ahb_spline[1](time),
    //                                    ahb_spline[2](time)};

    // Strahlkorper<Frame::Grid> horizon_a{2, 2, 1.0, std::move(aha_position)};
    // Strahlkorper<Frame::Grid> horizon_b{2, 2, 1.0, std::move(ahb_position)};

    // return std::make_pair(std::move(horizon_a), std::move(horizon_b));

    return control_system::TestHelpers::
        build_horizons_for_basic_control_systems<element_component>(
            time, runner, position_function, coord_map);
  };

  // Run the actual control system test.
  system_helper.run_control_system_test(runner, final_time, make_not_null(&gen),
                                        horizon_function, 2);
}
}  // namespace

// Currently the test takes a long time because the logic for the control system
// isn't the most optimal. This should be fixed in the near future.
// [[Timeout, 25]]
SPECTRE_TEST_CASE("Unit.ControlSystem.Systems.Compare",
                  "[ControlSystem][Unit]") {
  test_rotscaletrans_control_system<0, 3, 2>();
}
}  // namespace control_system
