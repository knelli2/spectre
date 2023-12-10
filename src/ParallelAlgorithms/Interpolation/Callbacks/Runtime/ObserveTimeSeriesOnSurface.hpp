// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <pup.h>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataBox/ValidateSelection.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Interpolation/Callbacks/Runtime/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace observers {
template <class Metavariables>
struct ObserverWriter;
}  // namespace observers
/// \endcond

namespace intrp {
namespace callbacks {
/// \cond
template <typename Target, typename TagsToObserve>
struct ObserveTimeSeriesOnSurface2;
/// \endcond

/// \brief post_interpolation_callback that outputs
/// a time series on a surface.
///
/// Uses:
/// - Metavariables
///   - `temporal_id`
/// - DataBox:
///   - `TagsToObserve`
///
/// Conforms to the intrp::protocols::PostInterpolationCallback protocol
///
/// For requirements on InterpolationTargetTag, see
/// intrp::protocols::InterpolationTargetTag
template <typename Target, typename... TagsToObserve>
struct ObserveTimeSeriesOnSurface2<Target, tmpl::list<TagsToObserve...>>
    : public Callback<Target> {
  static_assert(... and std::is_same_v<typename TagsToObserve::type, double>);

  /// \cond
  ObserveTimeSeriesOnSurface2() = default;
  explicit ObserveTimeSeriesOnSurface2(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(Callback, ObserveTimeSeriesOnSurface2);
  /// \endcond

  struct SubfileName {
    using type = std::string;
    static constexpr Options::String help = {
        "The name of the subfile inside the HDF5 file without an extension and "
        "without a preceding '/'."};
  };
  struct ValuesToObserve {
    using type = std::vector<std::string>;
    static constexpr Options::String help = {
        "List specifying the name of each value to observe."};
  };

  using options = tmpl::list<SubfileName, ValuesToObserve>;

  using available_tags_to_observe = tmpl::list<TagsToObserve...>;

  static constexpr Options::String help = {
      "Observe values (doubles) over a whole surface. These are usually "
      "reduced or integrated quantities."};

  ObserveTimeSeriesOnSurface2(std::string subfile_name,
                              const std::vector<std::string>& values_to_observe,
                              const Options::Context& context = {})
      : subfile_name_(std::move(subfile_name)) {
    db::validate_selection<available_tags_to_observe>(values_to_observe);
    for (const auto& value : values_to_observe) {
      values_to_observe_.insert(value);
    }
  }

  void pup(PUP::er& p) override {}

  template <typename Metavariables>
  static void apply(const db::Access& access,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const double time) {
    auto& proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);

    const std::unordered_set<std::string>& vars_to_observe =
        db::get<intrp::Tags::VarsToObserve>(access);

    std::vector<std::string> legend{};
    std::vector<double> data{};

    const size_t max_size = sizeof...(TagsToObserve) + 1;

    legend.reserve(max_size);
    data.reserve(max_size);

    legend.emplace_back("Time");
    legend.emplace_back(time);

    tmpl::for_each<available_tags_to_observe>(
        [&access, &vars_to_observe, &legend, &data](auto tag_v) {
          using Tag = tmpl::type_from<decltype(tag_v)>;
          const std::string tag_name = db::tag_name<Tag>();
          if (vars_to_observe.count(tag_name) == 1) {
            legend.emplace_back(tag_name);
            data.emplace_back(db::get<Tag>(access));
          }
        });

    // We call this on proxy[0] because the 0th element of a NodeGroup is
    // always guaranteed to be present.
    Parallel::threaded_action<
        observers::ThreadedActions::WriteReductionDataRow>(
        proxy[0], std::string{"/" + subfile_name_}, std::move(legend),
        std::make_tuple(std::move(data)));
  }

  const std::unordered_set<std::string>& observables() const override {
    return values_to_observe_;
  }

 private:
  std::string subfile_name_{};
  std::unordered_set<std::string> values_to_observe_{};
};

/// \cond
template <typename Target, typename... TagsToObserve>
PUP::able::PUP_ID ObserveTimeSeriesOnSurface2<
    Target, tmpl::list<TagsToObserve...>>::my_PUP_ID = 0;  // NOLINT
/// \endcond
}  // namespace callbacks
}  // namespace intrp
