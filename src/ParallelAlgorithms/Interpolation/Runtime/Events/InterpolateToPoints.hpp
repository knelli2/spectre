// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <pup.h>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/ObservationBox.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataBox/TagTraits.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/IdPair.hpp"
#include "Domain/ElementLogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "Options/Auto.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Callbacks/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Targets/Target.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/OptionalHelpers.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Events::Tags {
template <size_t Dim>
struct ObserverMesh;
}  // namespace Events::Tags
namespace intrp2 {
template <class Metavariables, bool IncludeDenseTriggers>
struct InterpolationTarget;
namespace Actions {
struct ReceiveVolumeTensors;
}  // namespace Actions
}  // namespace intrp2
/// \endcond

namespace intrp2::Events {
struct MarkAsInterpolation {
  // This will be used to initialize the Accesses that correspond to the
  // individual targets
  virtual std::unique_ptr<db::Access> initialize_target_element_box() const = 0;
};

struct ExampleOfTarget {
  static std::string name() { return "ExampleOfTarget"; }

  // The tag of the temporal ID to use.
  using temporal_id_tag = ::Tags::Time;

  // Whether or not you're allowed to create this target from options in the
  // EventsAnd{Dense}Triggers section.
  // QUESTION: Do we even want this? But isn't it necessary so we don't
  // accidentally create the control system events?
  static constexpr bool creatable_from_options = true;

  // The frame of this target. If it is `NoSuchType`, then the frame is chosen
  // at runtime. Otherwise it must be one of `Grid`, `Distorted`, `Inertial`.
  using frame = NoSuchType;
  // Struct containing a way to calculate the target points
  using points = Sphere;

  // List of callbacks that the target can possibly create at runtime. These
  // callbacks may or may not be called. It depends on if they are specific in
  // the input file
  using possible_runtime_callbacks = tmpl::list<>;
  // List of callbacks that are *always* called when the interpolation happens
  using compile_time_callbacks = tmpl::list<>;
};

template <typename TagsExclusivelyOnTarget, typename TagsFromVolume>
struct ExampleCallback {
  // Both these lists have simple and compute tags
  using tags_to_observe_on_target = tmpl::list<>;
  using non_observation_tags_on_target = tmpl::list<>;

  using volume_compute_tags = tmpl::list<>;
};

/*!
 * \brief Does an interpolation onto an InterpolationTargetTag by calling
 * Actions on the InterpolationTarget component.
 *
 * \note The `intrp2::TargetPoints::Sphere` target is handled specially because
 * it has the potential to be very slow due to it usually having the most points
 * out of all the stationary targets. An optimization for the future would be to
 * have each target be responsible for intelligently computing the
 * `block_logical_coordinates` for it's own points.
 */
template <typename Target, size_t Dim>
class InterpolateToPoints : public Event, MarkAsInterpolation {
  using temporal_id_tag = typename Target::temporal_id_tag;
  using TemporalId = typename temporal_id_tag::type;
  using PointsType = typename Target::points;
  using points_tags = typename PointsType::tags_on_target;
  using runtime_callbacks = typename Target::possible_runtime_callbacks;
  using compile_time_callbacks = typename Target::compile_time_callbacks;
  using all_callbacks = tmpl::append<runtime_callbacks, compile_time_callbacks>;

  using Frame = typename Target::frame;
  static constexpr bool use_runtime_frame = std::is_same_v<Frame, NoSuchType>;
  static constexpr bool use_runtime_callbacks =
      tmpl::size<runtime_callbacks>::value > 0;

  template <typename LocalCallback>
  struct get_target_tensors {
    using type = typename LocalCallback::tags_to_observe_on_target;
  };

  template <typename LocalCallback>
  struct get_non_observation_tags {
    using type = typename LocalCallback::non_observation_tags_on_target;
  };

  template <typename LocalCallback>
  struct get_volume_compute_tags {
    using type = typename LocalCallback::volume_compute_tags;
  };

  using all_target_tensors_to_observe = tmpl::remove_duplicates<tmpl::flatten<
      tmpl::transform<all_callbacks, get_target_tensors<tmpl::_1>>>>;

  using all_non_observation_tags_on_targets =
      tmpl::remove_duplicates<tmpl::flatten<
          tmpl::transform<all_callbacks, get_non_observation_tags<tmpl::_1>>>>;

  using common_tags_on_target =
      tmpl::list<intrp2::Tags::Frame, intrp2::Tags::Points<Dim>,
                 intrp2::Tags::VarsToObserve>;

  using all_tags_on_target =
      tmpl::remove_duplicates<tmpl::append<all_target_tensors_to_observe,
                                           all_non_observation_tags_on_targets,
                                           points_tags, common_tags_on_target>>;

 public:
  /// \cond
  explicit InterpolateToPoints(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(InterpolateToPoints);  // NOLINT
  /// \endcond

  static std::string name() { return Target::name(); }

  struct InputFrame {
    using type = std::string;
    static constexpr Options::String help = {
        "The frame of the points. This is one of 'Grid', 'Distorted', or "
        "'Inertial'. If you choose 'Distorted', make sure your domain has a "
        "distorted frame."};
  };

  struct Points {
    using type = PointsType;
    static constexpr Options::String help = {
        "Set of points to interpolate to."};
  };

 private:
  struct NaN {};

 public:
  struct InvalidPointFillValue {
    using type = Options::Auto<double, NaN>;
    static constexpr Options::String help = {
        "If a point is not in the computational domain, fill the variables "
        "that would have been interpolated to this point with this input "
        "value"};
  };

  struct Callbacks {
    static constexpr Options::String help = {"List of callbacks to run."};
    using type =
        std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>>;
    static size_t lower_bound_on_size() { return 1; }
  };

  using options = tmpl::flatten<
      tmpl::list<Points, InvalidPointFillValue,
                 tmpl::conditional_t<use_runtime_frame, tmpl::list<InputFrame>,
                                     tmpl::list<>>,
                 tmpl::conditional_t<use_runtime_callbacks,
                                     tmpl::list<Callbacks>, tmpl::list<>>>>;
  static constexpr Options::String help = {
      "Interpolate variables to a set of 'Points'."};

  InterpolateToPoints() = default;

  // Constructor for a compile-time frame. The frame won't be `NoSuchType` here
  // so it's ok to use `pretty_type::name` to get the string name of the frame.
  InterpolateToPoints(
      PointsType target, std::optional<double> invalid_fill_value,
      std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>>
          callbacks = std::vector<
              std::unique_ptr<intrp2::callbacks::Callback<Target>>>{},
      const Options::Context& context = {})
      : InterpolateToPoints(std::move(target), pretty_type::name<Frame>(),
                            std::move(callbacks), context) {}

  // Constructor for a runtime frame
  InterpolateToPoints(
      PointsType target, std::optional<double> invalid_fill_value,
      std::string frame,
      std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>>
          callbacks = std::vector<
              std::unique_ptr<intrp2::callbacks::Callback<Target>>>{},
      const Options::Context& context = {})
      : frame_(std::move(frame)),
        invalid_fill_value_(invalid_fill_value),
        target_(std::move(target)),
        callbacks_(std::move(callbacks)) {
    tensors_to_observe_.clear();
    for (const auto& callback : callbacks) {
      const std::unordered_set<std::string>& observables =
          callback->observables();
      tensors_to_observe_.insert(observables.begin(), observables.end());
    }

    // Compute the mapping between tensors to observe and the necessary volume
    // tensors to send
    compute_tensor_maps();

    for (const std::string& tensor_to_observe : tensors_to_observe_) {
      for (const std::string& volume_tag :
           observation_tags_to_volume_tags_.at(tensor_to_observe)) {
        volume_tensors_to_send_.insert(volume_tag);
      }
    }
  }

  using compute_tags_for_observation_box =
      tmpl::remove_duplicates<tmpl::flatten<
          tmpl::transform<all_callbacks, get_volume_compute_tags<tmpl::_1>>>>;

  using argument_tags = tmpl::list<::Tags::ObservationBox, temporal_id_tag,
                                   ::Events::Tags::ObserverMesh<Dim>>;

  template <typename DbTags, typename ComputeTagList,
            typename ParallelComponent, typename Metavariables>
  void operator()(const ObservationBox<DbTags, ComputeTagList>& box,
                  const TemporalId& temporal_id, const Mesh<Dim>& mesh,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ElementId<Dim>& array_index,
                  const ParallelComponent* const /*meta*/,
                  const ObservationValue& /*observation_value*/) const;

  using is_ready_argument_tags = tmpl::list<>;

  template <typename ArrayIndex, typename Component, typename Metavariables>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  // This will be used to initialize the Accesses that correspond to the
  // individual targets
  std::unique_ptr<db::Access> initialize_target_element_box() const override;

  // NOLINTNEXTLINE
  void pup(PUP::er& p);

  bool needs_evolved_variables() const override { return true; }

 private:
  std::string frame_{};
  std::optional<double> invalid_fill_value_{};
  PointsType target_{};
  std::vector<std::unique_ptr<intrp2::callbacks::Callback<Target>>>
      callbacks_{};
  std::unordered_set<std::string> tensors_to_observe_{};
  std::unordered_set<std::string> volume_tensors_to_send_{};

  template <typename ComputeTag, typename F>
  void for_compute_tag_arguments(const F&& f);

  void compute_tensor_maps();

  std::unordered_set<std::string> all_points_tags_{};
  std::unordered_map<std::string, std::vector<std::string>>
      non_obs_tags_to_volume_tags_{};
  std::unordered_map<std::string, std::vector<std::string>>
      observation_tags_to_volume_tags_{};
};

template <typename Target, size_t Dim>
template <typename ComputeTag, typename F>
void InterpolateToPoints<Target, Dim>::for_compute_tag_arguments(const F&& f) {
  using argument_tags = typename ComputeTag::argument_tags;

  tmpl::for_each<argument_tags>([&f](auto tag_v) { f(tag_v); });
}

template <typename Target, size_t Dim>
void InterpolateToPoints<Target, Dim>::compute_tensor_maps() {
  // Create list of tags that are specific to the points we are using and have
  // no relation to the volume quantities whatsoever. These tags cannot be
  // observed.
  tmpl::for_each<points_tags>([this](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;

    all_points_tags_.insert(db::tag_name<Tag>());

    if constexpr (db::is_compute_tag_v<Tag>) {
      for_compute_tag_arguments<Tag>([this](auto argument_tag_v) {
        using ArgumentTag = tmpl::type_from<decltype(argument_tag_v)>;

        all_points_tags_.insert(db::tag_name<Tag>());
      });
    }
  });

  // Create map from the non-observable tags to the required volume tags
  // necessary. This is needed in case an actual observable tag depends on one
  // of these non-observable tags, which in turn needs volume vars. However, we
  // only specify observable tags in the input file.
  tmpl::for_each<all_non_observation_tags_on_targets>([this](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;

    const std::string tag_name = db::tag_name<Tag>();

    // If it's a simple tag, then we assume it's a volume tag
    if constexpr (db::is_simple_tag_v<Tag>) {
      non_obs_tags_to_volume_tags_[tag_name] =
          std::vector<std::string>{tag_name};
    } else {
      static_assert(db::is_compute_tag_v<Tag>);

      non_obs_tags_to_volume_tags_[tag_name];
      non_obs_tags_to_volume_tags_.at(tag_name).reserve(
          tmpl::size<typename Tag::argument_tags>::value);

      // Go through all arguments of the compute tag. If an argument tag is in
      // the points tags, then don't add it to the map.
      for_compute_tag_arguments<Tag>([this](auto argument_tag_v) {
        using ArgumentTag = tmpl::type_from<decltype(argument_tag_v)>;

        std::string argument_tag_name = db::tag_name<ArgumentTag>();

        if (all_points_tags_.count(argument_tag_name) == 0) {
          non_obs_tags_to_volume_tags_.at(tag_name).emplace_back(
              std::move(argument_tag_name));
        }
      });

      // QUESTION: Shouldn't happen?
      if (non_obs_tags_to_volume_tags_.at(tag_name).size() == 0) {
        non_obs_tags_to_volume_tags_.erase(tag_name);
      }
    }
  });

  // Create the actual map that goes from the tags users can specify in the
  // input file to volume variables that need to be interpolated to the target.
  tmpl::for_each<all_target_tensors_to_observe>([this](auto tag_v) {
    using Tag = tmpl::type_from<decltype(tag_v)>;

    const std::string tag_name = db::tag_name<Tag>();

    if constexpr (db::is_simple_tag_v<Tag>) {
      observation_tags_to_volume_tags_[tag_name] =
          std::vector<std::string>{tag_name};
    } else {
      static_assert(db::is_compute_tag_v<Tag>);

      observation_tags_to_volume_tags_[tag_name];

      // Go through all arguments of the compute tag. If an argument tag is in
      // the points tags, then don't add it to the map.
      for_compute_tag_arguments<Tag>([this](auto argument_tag_v) {
        using ArgumentTag = tmpl::type_from<decltype(argument_tag_v)>;

        std::string argument_tag_name = db::tag_name<ArgumentTag>();

        // Don't add tags that are in the points tags
        if (all_points_tags_.count(argument_tag_name) == 1) {
          return;
        }

        // Check if this argument was from the non-observable tags. If it was,
        // then add all of those volume tags to the map for this tag
        if (non_obs_tags_to_volume_tags_.count(argument_tag_name) == 1) {
          observation_tags_to_volume_tags_.at(tag_name).insert(
              observation_tags_to_volume_tags_.at(tag_name).end(),
              non_obs_tags_to_volume_tags_.at(argument_tag_name).begin(),
              non_obs_tags_to_volume_tags_.at(argument_tag_name).end());
        } else {
          // Just a simple tag, so this is a volume tag
          observation_tags_to_volume_tags_.at(tag_name).emplace_back(
              std::move(argument_tag_name));
        }
      });
    }
  });
}

template <typename Target, size_t Dim>
template <typename DbTags, typename ComputeTagList, typename ParallelComponent,
          typename Metavariables>
void InterpolateToPoints<Target, Dim>::operator()(
    const ObservationBox<DbTags, ComputeTagList>& box,
    const TemporalId& temporal_id, const Mesh<Dim>& mesh,
    Parallel::GlobalCache<Metavariables>& cache,
    const ElementId<Dim>& array_index, const ParallelComponent* const /*meta*/,
    const ObservationValue& /*observation_value*/) const {
  // TODO: Fix the time argument here. Need something like
  // InterpolationTarget_detail::get_temporal_id_value
  const std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>
      block_logical_coords = target_->block_logical_coordinates(
          box, cache, temporal_id.time(), frame_);

  const std::vector<ElementId<Dim>> element_ids{{array_index}};
  const auto element_coord_holders =
      element_logical_coordinates(element_ids, block_logical_coords);

  if (element_coord_holders.count(array_index) == 0) {
    // There are no target points in this element, so we don't need
    // to do anything.
    return;
  }

  // There are points in this element, so interpolate to them and
  // send the interpolated data to the target.
  const auto& element_coord_holder = element_coord_holders.at(array_index);

  // string is the name of the tensor, std::vector is for each component of the
  // tensor, and the DataVector is the actual data
  std::unordered_map<std::string, std::vector<DataVector>> interpolated_vars{};

  // 2. Set up interpolator
  intrp::Irregular<Dim> interpolator(
      mesh, element_coord_holder.element_logical_coords);
  const size_t number_of_points_in_this_element =
      element_coord_holder.element_logical_coords.get(0).size();

  using all_obs_box_tensors =
      typename ObservationBox<DbTags, ComputeTagList>::tags_list;

  // Go through and interpolate all necessary volume tensor components
  tmpl::for_each<all_obs_box_tensors>(
      [this, &box, &interpolated_vars, &interpolator,
       &number_of_points_in_this_element](auto tag_v) {
        using Tag = tmpl::type_from<decltype(tag_v)>;

        const std::string& tag_name = db::tag_name<Tag>();

        // If we don't need to send this tensor, then skip it
        if (not volume_tensors_to_send_.count(tag_name) == 1) {
          return;
        }

        const auto& tensor = db::get<Tag>(box);

        const size_t num_tensor_indices = Tag::type::num_tensor_indices;
        interpolated_vars[tag_name] = std::vector<DataVector>{
            num_tensor_indices, DataVector{number_of_points_in_this_element}};

        for (size_t i = 0; i < num_tensor_indices; i++) {
          interpolator.interpolate(
              make_not_null(&interpolated_vars.at(tag_name)[i]), tensor.get(i));
        }
      });

  auto& receiver_proxy =
      Parallel::get_parallel_component<InterpolationTarget<Metavariables>>(
          cache)[name()];
  Parallel::simple_action<intrp2::Actions::ReceiveVolumeTensors>(
      receiver_proxy, temporal_id, std::move(interpolated_vars),
      element_coord_holder.offsets);
}

template <typename Target, size_t Dim>
std::unique_ptr<db::Access>
InterpolateToPoints<Target, Dim>::initialize_target_element_box() const {
  using BoxType = db::compute_databox_type<all_tags_on_target>;

  BoxType typed_box{};

  db::mutate_apply<common_tags_on_target, tmpl::list<>>(
      [this](
          const gsl::not_null<std::string*> frame,
          const gsl::not_null<tnsr::I<DataVector, Dim, Frame::NoFrame>*> points,
          const gsl::not_null<std::unordered_set<std::string>*> vars_to_observe,
          const gsl::not_null<double*> invalid_points_fill_value) {
        *frame = frame_;
        *points = target_->target_points_no_frame();
        *vars_to_observe = tensors_to_observe_;
        *invalid_points_fill_value = invalid_fill_value_.value_or(
            std::numeric_limits<double>::quiet_NaN());
      },
      make_not_null(&typed_box));

  return std::make_unique<BoxType>(std::move(typed_box));
}

template <typename Target, size_t Dim>
void InterpolateToPoints<Target, Dim>::pup(PUP::er& p) {
  p | frame_;
  p | invalid_fill_value_;
  p | target_;
  p | callbacks_;
  p | tensors_to_observe_;
}

/// \cond
template <typename Target, size_t Dim>
PUP::able::PUP_ID InterpolateToPoints<Target, Dim>::my_PUP_ID = 0;  // NOLINT
/// \endcond

}  // namespace intrp2::Events
