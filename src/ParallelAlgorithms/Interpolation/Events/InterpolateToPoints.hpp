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
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/IdPair.hpp"
#include "Domain/ElementLogicalCoordinates.hpp"
#include "Domain/Structure/BlockId.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/Interpolation/IrregularInterpolant.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Options/String.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Targets/RuntimeTargets/Target.hpp"
#include "Utilities/OptionalHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Events::Tags {
template <size_t Dim>
struct ObserverMesh;
}  // namespace Events::Tags
/// \endcond

namespace intrp::Events {
/// \cond
template <size_t Dim, typename SourceVarTags>
class InterpolateToPoints;
/// \endcond

/*!
 * \brief Does an interpolation onto an InterpolationTargetTag by calling
 * Actions on the InterpolationTarget component.
 *
 * \note The `intrp::TargetPoints::Sphere` target is handled specially because
 * it has the potential to be very slow due to it usually having the most points
 * out of all the stationary targets. An optimization for the future would be to
 * have each target be responsible for intelligently computing the
 * `block_logical_coordinates` for it's own points.
 */
template <size_t Dim, typename... SourceVarTags>
class InterpolateToPoints<Dim, tmpl::list<SourceVarTags...>> : public Event {
 public:
  /// \cond
  explicit InterpolateToPoints(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(InterpolateToPoints);  // NOLINT
  /// \endcond

  struct InputFrame {
    using type = std::string;
    static constexpr Options::String help = {
        "The frame of the points. This is one of 'Grid', 'Distorted', or "
        "'Inertial'. If you choose 'Distorted', make sure your domain has a "
        "distorted frame."};
  };

  struct Points {
    using type = std::unique_ptr<intrp::Targets::Target<Dim>>;
    static constexpr Options::String help = {
        "Set of points to interpolate to."};
  };

  struct VariablesToObserve {
    static constexpr Options::String help = "Subset of variables to observe";
    using type = std::vector<std::string>;
    static size_t lower_bound_on_size() { return 1; }
  };

  using options = tmpl::list<InputFrame, Points VariablesToObserve>;
  static constexpr Options::String help = {
      "Interpolate to a set of 'Points' in a specific 'Frame'."};

  InterpolateToPoints() = default;

  InterpolateToPoints(std::string frame,
                      std::unique_ptr<intrp::Targets::Target<Dim>> target,
                      const std::vector<std::string>& variables_to_observe,
                      const Options::Context& context = {})
      : frame_(frame),
        target_(std::move(target)),
        variables_to_observe_([&context, &variables_to_observe]() {
          std::unordered_set<std::string> result{};
          for (const auto& tensor : variables_to_observe) {
            if (result.count(tensor) != 0) {
              PARSE_ERROR(
                  context,
                  "Listed variable '"
                      << tensor
                      << "' more than once in list of variables to observe.");
            }
            result.insert(tensor);
          }
          return result;
        }()) {
    std::unordered_set<std::string> valid_tensors{};
    tmpl::for_each<available_tags_to_observe>([&valid_tensors](auto tag_v) {
      using tag = tmpl::type_from<decltype(tag_v)>;
      valid_tensors.insert(detail::name<tag>());
    });

    for (const auto& name : variables_to_observe_) {
      if (valid_tensors.count(name) != 1) {
        PARSE_ERROR(
            context,
            name << " is not an available variable. Available variables:\n"
                 << valid_tensors);
      }
    }
  }

  using compute_tags_for_observation_box = tmpl::list<>;

  using argument_tags =
      tmpl::list<::Tags::DataBox, Tags::Time, ::Events::Tags::ObserverMesh<Dim>,
                 SourceVarTags...>;

  template <typename DbTags, typename ParallelComponent, typename Metavariables>
  void operator()(const db::DataBox<DbTags>& box, const double time,
                  const Mesh<Dim>& mesh,
                  const typename SourceVarTags::type&... source_vars_input,
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

  // NOLINTNEXTLINE
  void pup(PUP::er& p);

  bool needs_evolved_variables() const override { return true; }

 private:
  std::string frame_{};
  std::unique_ptr<intrp::Targets::Target<Dim>> target_{};
  std::unordered_set<std::string> variables_to_observe_{};
};

template <size_t Dim, typename... SourceVarTags>
template <typename DbTags, typename ParallelComponent, typename Metavariables>
void InterpolateToPoints<Dim, tmpl::list<SourceVarTags...>>::operator()(
    const db::DataBox<DbTags>& box, const double time, const Mesh<Dim>& mesh,
    const typename SourceVarTags::type&... source_vars_input,
    Parallel::GlobalCache<Metavariables>& cache,
    const ElementId<Dim>& array_index, const ParallelComponent* const /*meta*/,
    const ObservationValue& /*observation_value*/) const {
  const std::vector<std::optional<
      IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>
      block_logical_coords =
          target_->block_logical_coordinates(box, cache, time, frame_);

  const std::vector<ElementId<Dim>> element_ids{{array_index}};
  const auto element_coord_holders =
      element_logical_coordinates(element_ids, block_logical_coords);

  if (element_coord_holders.count(array_index) == 0) {
    // There are no target points in this element, so we don't need
    // to do anything.
    return;
  }

  std::vector<double> interpolated_vars{};
  // This is larger than we need if we are only observing some tensors, but
  // that's not a big deal and calculating the correct size is nontrivial.
  interpolated_vars.reserve(alg::accumulate(
      std::initializer_list<size_t>{
          mesh.number_of_grid_points() *
          std::decay_t<decltype(
              value(typename SourceVarTags::type{}))>::size()...},
      0_st));

  // There are points in this element, so interpolate to them and
  // send the interpolated data to the target.  This is done
  // in several steps:
  const auto& element_coord_holder = element_coord_holders.at(array_index);

  // 2. Set up interpolator
  intrp::Irregular<Dim> interpolator(
      mesh, element_coord_holder.element_logical_coords);

  // 3. Interpolate and send interpolated data to target
  // auto& receiver_proxy = Parallel::get_parallel_component<
  //     InterpolationTarget<Metavariables, InterpolationTargetTag>>(cache);
  // Parallel::simple_action<
  //     Actions::InterpolationTargetVarsFromElement<InterpolationTargetTag>>(
  //     receiver_proxy,
  //     std::vector<Variables<
  //         typename InterpolationTargetTag::vars_to_interpolate_to_target>>(
  //         {interpolator.interpolate(interp_vars)}),
  //     block_logical_coords,
  //     std::vector<std::vector<size_t>>({element_coord_holder.offsets}),
  //     time);
}

template <size_t Dim, typename... SourceVarTags>
void InterpolateToPoints<Dim, tmpl::list<SourceVarTags...>>::pup(PUP::er& p) {
  p | frame_;
  p | target_;
  p | variables_to_observe_;
}

/// \cond
template <size_t Dim, typename... SourceVarTags>
PUP::able::PUP_ID InterpolateToPoints<Dim,
                                      tmpl::list<SourceVarTags...>>::my_PUP_ID =
    0;  // NOLINT
/// \endcond

}  // namespace intrp::Events
