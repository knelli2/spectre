// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <random>
#include <string>

#include "ApparentHorizons/ObserveCenters.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Matrix.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/Observer/Initialize.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Tags.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace {

// This is needed to template ObserveCenters, but it doesn't have to be an
// actual InterpolationTargetTag. This is only used for the name of the subfile
// written
struct AhA {};

template <typename Metavariables>
struct MockObserverWriter {
  using component_being_mocked = observers::ObserverWriter<Metavariables>;
  using replace_these_simple_actions = tmpl::list<>;
  using with_these_simple_actions = tmpl::list<>;

  using const_global_cache_tags =
      tmpl::list<observers::Tags::ReductionFileName>;

  using initialize_action_list =
      tmpl::list<::Actions::SetupDataBox,
                 observers::Actions::InitializeWriter<Metavariables>>;
  using initialization_tags =
      Parallel::get_initialization_tags<initialize_action_list>;

  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockNodeGroupChare;
  using array_index = int;

  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        initialize_action_list>>;
};

struct TestMetavars {
  using observed_reduction_data_tags = tmpl::list<>;

  enum class Phase { Initialization, WriteData, Exit };

  using component_list = tmpl::list<MockObserverWriter<TestMetavars>>;
};

using FoTPtr = std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>;

SPECTRE_TEST_CASE("Unit.ApparentHorizons.ObserveCenters",
                  "[Unit][ApparentHorizons]") {
  using observer = MockObserverWriter<TestMetavars>;
  MAKE_GENERATOR(gen);
  std::uniform_real_distribution<double> center_dist{-10.0, 10.0};
  const std::string filename = "ObserveCentersTest_Output";

  // clean up just in case
  if (file_system::check_if_file_exists(filename + ".h5")) {
    file_system::rm(filename + ".h5", true);
  }

  // set up runner and stuff
  ActionTesting::MockRuntimeSystem<TestMetavars> runner{{filename}};
  runner.set_phase(TestMetavars::Phase::Initialization);
  ActionTesting::emplace_nodegroup_component<observer>(make_not_null(&runner));
  auto& cache = ActionTesting::cache<observer>(runner, 0);

  // the initialization actions
  for (size_t i = 0; i < 2; ++i) {
    ActionTesting::next_action<observer>(make_not_null(&runner), 0);
  }

  runner.set_phase(TestMetavars::Phase::WriteData);

  const auto make_center = [&gen, &center_dist]() -> std::array<double, 3> {
    return make_with_random_values<std::array<double, 3>>(
        make_not_null(&gen), center_dist, std::array<double, 3>{});
  };

  // Lists of centers so we can check that the correct centers were written
  std::vector<std::array<double, 3>> grid_centers{};
  std::vector<std::array<double, 3>> inertial_centers{};

  db::DataBox<tmpl::list<StrahlkorperTags::Strahlkorper<Frame::Grid>,
                         StrahlkorperTags::Strahlkorper<Frame::Inertial>>>
      box{};

  const auto update_stored_centers = [&make_center, &grid_centers,
                                      &inertial_centers, &box]() {
    const auto grid_center = make_center();
    const auto inertial_center = make_center();

    grid_centers.push_back(grid_center);
    inertial_centers.push_back(inertial_center);

    db::mutate<StrahlkorperTags::Strahlkorper<Frame::Grid>,
               StrahlkorperTags::Strahlkorper<Frame::Inertial>>(
        make_not_null(&box),
        [&grid_center, &inertial_center](
            gsl::not_null<Strahlkorper<Frame::Grid>*> box_grid_horizon,
            gsl::not_null<Strahlkorper<Frame::Inertial>*>
                box_inertial_horizon) {
          // We only care about the centers of the strahlkorper so the
          // ell and radius are arbitrary
          *box_grid_horizon = Strahlkorper<Frame::Grid>{2, 1.0, grid_center};
          *box_inertial_horizon =
              Strahlkorper<Frame::Inertial>{2, 1.0, inertial_center};
        });
  };

  // times to write
  const std::vector<double> times{0.0, 0.1, 0.2, 0.3, 0.4, 0.5};

  // write some data
  for (size_t i = 0; i < times.size(); i++) {
    update_stored_centers();

    ah::callbacks::ObserveCenters<AhA>::apply(box, cache, times[i]);

    size_t num_threaded_actions =
        ActionTesting::number_of_queued_threaded_actions<observer>(runner, 0);
    CHECK(num_threaded_actions == 1);
    ActionTesting::invoke_queued_threaded_action<observer>(
        make_not_null(&runner), 0);
  }

  // scoped to close file
  {
    h5::H5File<h5::AccessType::ReadOnly> read_file{filename + ".h5"};
    // These must be the same as in ObserveCenters
    const std::vector<std::string> compare_legend{
        {"Time", "GridCenter_X", "GridCenter_Y", "GridCenter_Z",
         "InertialCenter_X", "InertialCenter_Y", "InertialCenter_Z"}};
    const std::string subfile_name =
        "/ApparentHorizons/" + pretty_type::name<AhA>() + "_Centers";

    const auto& dataset = read_file.get<h5::Dat>(subfile_name);
    const Matrix data = dataset.get_data();
    const std::vector<std::string>& legend = dataset.get_legend();

    // Check legend is correct
    for (size_t i = 0; i < legend.size(); i++) {
      CHECK(legend[i] == compare_legend[i]);
    }

    // Check proper number of times were written
    CHECK(data.rows() == times.size());

    // Check centers
    for (size_t i = 0; i < times.size(); i++) {
      CHECK(data(i, 0) == times[i]);

      const std::array<double, 3>& grid_center = grid_centers[i];
      const std::array<double, 3>& inertial_center = inertial_centers[i];
      for (size_t j = 0; j < grid_center.size(); j++) {
        // Grid center is columns 2-4
        CHECK(data(i, j + 1) == gsl::at(grid_center, j));
        // Inertial center is columns 5-7
        CHECK(data(i, j + 4) == gsl::at(inertial_center, j));
      }
    }
  }

  // clean up now that we are done
  if (file_system::check_if_file_exists(filename + ".h5")) {
    file_system::rm(filename + ".h5", true);
  }
}
}  // namespace
