// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Initialization/Tags.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/Tags.hpp"
#include "IO/Logging/Tags.hpp"
#include "IO/Logging/Verbosity.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/LtsTimeStepper.hpp"
#include "Utilities/Rational.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Cce {
namespace Actions {

template <typename Metavariables>
struct InitializeCcmNextTime {
  using simple_tags_from_options = tmpl::list<>;
  using const_global_cache_tags =
      tmpl::list<logging::Tags::Verbosity<Cce::OptionTags::Cce>>;
  using mutable_global_cache_tags = tmpl::list<>;
  using return_tags = tmpl::list<Cce::Tags::NextCcmTime>;
  using argument_tags = tmpl::list<::Tags::TimeStepId>;
  using compute_tags = tmpl::list<>;
  using simple_tags =
      tmpl::append<return_tags, typename Metavariables::ccm_psi0>;

  static void apply(gsl::not_null<TimeStepId*> next_ccm_time,
                    const TimeStepId& current_time) {
    *next_ccm_time = current_time;
  }
};

}  // namespace Actions
}  // namespace Cce
