// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <random>
#include <tuple>
#include <utility>

#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Helpers/PointwiseFunctions/PostNewtonian/BinaryTrajectories.hpp"

namespace {
// must match value in BinaryTrajectories.py
constexpr double initial_separation{15.366};
const std::array<double, 3> initial_velocity{0.1, -0.2, 0.3};

std::array<double, 3> positions1(const double time) {
  const BinaryTrajectories expected{initial_separation, initial_velocity};
  return expected.positions(time).first;
}

std::array<double, 3> positions2(const double time) {
  const BinaryTrajectories expected{initial_separation, initial_velocity};
  return expected.positions(time).second;
}
}  // namespace

SPECTRE_TEST_CASE("Test.TestHelpers.PostNewtonian.BinaryTrajectories",
                  "[PointwiseFunctions][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env(
      "Helpers/Tests/PointwiseFunctions/PostNewtonian/");
  const BinaryTrajectories expected{initial_separation, initial_velocity};

  pypp::check_with_random_values<1>(
      &BinaryTrajectories::separation,
      BinaryTrajectories{initial_separation, initial_velocity},
      "BinaryTrajectories", {"separation"}, {{{-100.0, 100.0}}},
      std::make_tuple(initial_separation), initial_separation);
  pypp::check_with_random_values<1>(
      &BinaryTrajectories::orbital_frequency,
      BinaryTrajectories{initial_separation, initial_velocity},
      "BinaryTrajectories", {"orbital_frequency"}, {{{-100.0, 100.0}}},
      std::make_tuple(initial_separation), initial_separation);
  pypp::check_with_random_values<1>(&positions1, "BinaryTrajectories",
                                    "positions1", {{{-100.0, 100.0}}},
                                    initial_separation);
  pypp::check_with_random_values<1>(&positions2, "BinaryTrajectories",
                                    "positions2", {{{-100.0, 100.0}}},
                                    initial_separation);
}
