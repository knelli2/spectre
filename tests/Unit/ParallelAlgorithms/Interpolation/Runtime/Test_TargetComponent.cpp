// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>

#include "Options/String.hpp"
#include "Parallel/Phase.hpp"
#include "Utilties/TMPL.hpp"

namespace {
// blah

struct TestMetavariables {
  using component_list = tmpl::list<>;

  static constexpr Options::String help{
      "An executable for testing the array allocation of the "
      "InterpolationTarget parallel component."};

  static constexpr std::array<Parallel::Phase, 3> default_phase_order{
      {Parallel::Phase::Initialization, Parallel::Phase::Testing,
       Parallel::Phase::Exit}};
}
}  // namespace
