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
#include "DataStructures/DataBox/TagTraits.hpp"
#include "DataStructures/DataBox/ValidateSelection.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "Options/Context.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Metafunctions.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Tags.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace observers {
template <class Metavariables>
struct ObserverWriter;
}  // namespace observers
/// \endcond

namespace intrp2 {
namespace callbacks {
/// \cond
template <typename Target, typename TagsToObserve, typename NonObservationTags,
          typename VolumeComputeTags>
struct ObserveTimeSeriesOnSurface;
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
/// Conforms to the intrp2::protocols::PostInterpolationCallback protocol
///
/// For requirements on InterpolationTargetTag, see
/// intrp2::protocols::InterpolationTargetTag
// TODO: Rename to ObserveTimeSeriesOnTarget
template <typename Target, typename... TagsToObserve,
          typename NonObservationTags, typename VolumeComputeTags>
struct ObserveTimeSeriesOnSurface<Target, tmpl::list<TagsToObserve...>,
                                  NonObservationTags, VolumeComputeTags>
    : public Callback<Target>, public tt::ConformsTo<protocols::Callback> {
  static_assert(tmpl2::flat_all_v<
                std::is_same_v<typename TagsToObserve::type, double>...>);

  using tags_to_observe_on_target = tmpl::list<TagsToObserve...>;
  using non_observation_tags_on_target = NonObservationTags;
  using volume_compute_tags = VolumeComputeTags;

 private:
  using simple_tags_to_observe =
      metafunctions::simple_tags_from_mixed_tags<tags_to_observe_on_target>;

 public:
  /// \cond
  ObserveTimeSeriesOnSurface() = default;
  explicit ObserveTimeSeriesOnSurface(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  // NOLINTNEXTLINE
  WRAPPED_PUPable_decl_base_template(Callback<Target>,
                                     ObserveTimeSeriesOnSurface);
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

  static constexpr Options::String help = {
      "Observe values (doubles) over a whole surface. These are usually "
      "reduced or integrated quantities."};

  ObserveTimeSeriesOnSurface(std::string subfile_name,
                             const std::vector<std::string>& values_to_observe,
                             const Options::Context& context = {})
      : subfile_name_(std::move(subfile_name)) {
    db::validate_selection<tags_to_observe_on_target>(values_to_observe,
                                                      context);
    for (const auto& value : values_to_observe) {
      values_to_observe_.insert(value);
    }
  }

  void pup(PUP::er& p) override {
    p | subfile_name_;
    p | values_to_observe_;
  }

  template <typename Metavariables>
  void apply(const db::Access& access,
             Parallel::GlobalCache<Metavariables>& cache, const double time) {
    auto& proxy = Parallel::get_parallel_component<
        observers::ObserverWriter<Metavariables>>(cache);

    std::vector<std::string> legend{};
    std::vector<double> data{};

    const size_t max_size = sizeof...(TagsToObserve) + 1;

    legend.reserve(max_size);
    data.reserve(max_size);

    legend.emplace_back("Time");
    data.emplace_back(time);

    tmpl::for_each<simple_tags_to_observe>(
        [this, &access, &legend, &data](auto tag_v) {
          using Tag = tmpl::type_from<decltype(tag_v)>;
          const std::string tag_name = db::tag_name<Tag>();
          if (values_to_observe_.contains(tag_name)) {
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
// NOLINTBEGIN
template <typename Target, typename... TagsToObserve,
          typename NonObservationTags, typename VolumeComputeTags>
PUP::able::PUP_ID ObserveTimeSeriesOnSurface<
    Target, tmpl::list<TagsToObserve...>, NonObservationTags,
    VolumeComputeTags>::my_PUP_ID = 0;
// NOLINTEND
/// \endcond
}  // namespace callbacks
}  // namespace intrp2
