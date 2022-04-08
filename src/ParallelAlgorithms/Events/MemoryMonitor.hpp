// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <optional>
#include <pup.h>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "IO/Observer/GetSectionObservationKey.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/Tags.hpp"
#include "IO/Observer/TypeOfObservation.hpp"
#include "Options/Options.hpp"
#include "Parallel/ArrayIndex.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/Serialize.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Numeric.hpp"
#include "Utilities/OptionalHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Parallel::Algorithms {
struct Array;
struct Group;
struct Nodegroup;
struct Singleton;
}  // namespace Parallel::Algorithms
/// \endcond

namespace Events {

namespace Tags {
template <typename Component>
struct ArraySectionIdTag : db::SimpleTag {
  using type = std::string;
  static std::string name() {
    return std::is_same_v<Component, void> ? "MutableGlobalCache"
                                           : pretty_type::name<Component>();
  }
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<>;
  static std::string create_from_options() { return name(); }
};
struct DesignatedNode : db::SimpleTag {
  using type = bool;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<>;
  static bool create_from_options() { return false; }
};
struct DesignatedProc : db::SimpleTag {
  using type = bool;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<>;
  static bool create_from_options() { return false; }
};
}  // namespace Tags

namespace Actions {
namespace detail {
template <typename Component>
using array_section_id_tag = ::Events::Tags::ArraySectionIdTag<Component>;
}  // namespace detail

template <typename Metavariables>
struct InitializeMemoryMonitor {
  using section_ids =
      tmpl::transform<typename Metavariables::component_list,
                      tmpl::bind<detail::array_section_id_tag, tmpl::_1>>;

  using observer_tags =
      tmpl::transform<section_ids,
                      tmpl::bind<observers::Tags::ObservationKey, tmpl::_1>>;

  using initialization_tags =
      tmpl::flatten<tmpl::list<section_ids, Events::Tags::DesignatedNode,
                               Events::Tags::DesignatedProc, observer_tags>>;

  using initialization_tags_to_keep = initialization_tags;

  using simple_tags = tmpl::list<>;

  using compute_tags = tmpl::list<>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    tmpl::for_each<section_ids>([&box](auto tag_v) {
      using section_id = tmpl::type_from<decltype(tag_v)>;
      using tag = observers::Tags::ObservationKey<section_id>;
      std::optional<std::string> name{};
      name = section_id::name();
      ::Initialization::mutate_assign<tmpl::list<tag>>(make_not_null(&box),
                                                       std::move(name));
    });

    return std::make_tuple(std::move(box));
  }
};
}  // namespace Actions

template <typename Proxy, typename TupleData>
void singleton_memory_monitor(Proxy& proxy, const std::string& filename,
                              const std::vector<std::string>& legend,
                              const TupleData& data) {
  Parallel::threaded_action<observers::ThreadedActions::WriteReductionDataRow>(
      proxy[0], filename, legend, data);
}

template <typename ObservationValueTag, typename Metavariables, size_t Index>
class MemoryMonitor : public Event {
 private:
  using ReductionData = tmpl::wrap<
      tmpl::list<
          // Time
          Parallel::ReductionDatum<double, funcl::AssertEqual<>>,
          // Node/Proc
          Parallel::ReductionDatum<int, funcl::AssertEqual<>>,
          // Number of elements
          Parallel::ReductionDatum<size_t, funcl::Plus<>>,
          // Memory usage in MB
          Parallel::ReductionDatum<double, funcl::Plus<>, funcl::Divides<>,
                                   std::index_sequence<2>>>,
      Parallel::ReductionData>;

 public:
  explicit MemoryMonitor(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(MemoryMonitor);  // NOLINT

  using options = tmpl::list<>;

  static std::string name() {
    using component_list = typename Metavariables::component_list;
    static_assert(Index <= tmpl::size<component_list>::value);
    using ComponentName =
        tmpl::conditional_t<tmpl::size<component_list>::value == Index, void,
                            tmpl::at<component_list, tmpl::size_t<Index>>>;
    return "MemoryMonitor" + pretty_type::name<ComponentName>();
  }

  static constexpr Options::String help =
      "Observe memory usage of each Parallel Component";

  MemoryMonitor() = default;

  using observed_reduction_data_tags =
      observers::make_reduction_data_tags<tmpl::list<ReductionData>>;

  using compute_tags_for_observation_box = tmpl::list<>;

  using argument_tags = tmpl::list<ObservationValueTag, ::Tags::DataBox>;

  template <typename DbTags, typename ArrayIndex, typename ParallelComponent>
  void operator()(const typename ObservationValueTag::type& observation_value,
                  const db::DataBox<DbTags>& box,
                  Parallel::GlobalCache<Metavariables>& cache,
                  const ArrayIndex& array_index,
                  const ParallelComponent* const /*meta*/) const;

  using observation_registration_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTagsList>
  std::optional<
      std::pair<observers::TypeOfObservation, observers::ObservationKey>>
  get_observation_type_and_key_for_registration(
      const db::DataBox<DbTagsList>& box) const {
    using component_list = typename Metavariables::component_list;
    static_assert(Index <= tmpl::size<component_list>::value);
    using ComponentName =
        tmpl::conditional_t<tmpl::size<component_list>::value == Index, void,
                            tmpl::at<component_list, tmpl::size_t<Index>>>;

    const std::optional<std::string> section_observation_key =
        observers::get_section_observation_key<
            Tags::ArraySectionIdTag<ComponentName>>(box);
    if (not section_observation_key.has_value()) {
      return std::nullopt;
    }
    return {{observers::TypeOfObservation::Reduction,
             observers::ObservationKey(
                 subfile_path_ + section_observation_key.value() + ".dat")}};
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return true; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

 private:
  std::string subfile_path_{"/MemoryMonitors/"};
};

/// \cond
template <typename ObservationValueTag, typename Metavariables, size_t Index>
MemoryMonitor<ObservationValueTag, Metavariables, Index>::MemoryMonitor(
    CkMigrateMessage* msg)
    : Event(msg) {}

template <typename ObservationValueTag, typename Metavariables, size_t Index>
template <typename DbTags, typename ArrayIndex, typename ParallelComponent>
void MemoryMonitor<ObservationValueTag, Metavariables, Index>::operator()(
    const typename ObservationValueTag::type& observation_value,
    const db::DataBox<DbTags>& box, Parallel::GlobalCache<Metavariables>& cache,
    const ArrayIndex& array_index,
    const ParallelComponent* const /*meta*/) const {
  using component_list = typename Metavariables::component_list;
  static_assert(Index <= tmpl::size<component_list>::value);
  using ComponentName =
      tmpl::conditional_t<tmpl::size<component_list>::value == Index, void,
                          tmpl::at<component_list, tmpl::size_t<Index>>>;

  // Skip observation on elements that are not part of a section
  const std::optional<std::string> section_observation_key =
      observers::get_section_observation_key<
          Tags::ArraySectionIdTag<ComponentName>>(box);
  if (not section_observation_key.has_value()) {
    return;
  }

  std::vector<std::string> legend{
      db::tag_name<ObservationValueTag>(), "Node/Proc",
      "Number of elements on Node/Proc", "Memory Usage (MB)"};

  const bool designated_proc = db::get<Tags::DesignatedProc>(box);
  const bool designated_node = db::get<Tags::DesignatedNode>(box);

  auto& local_observer = *Parallel::local_branch(
      Parallel::get_parallel_component<observers::Observer<Metavariables>>(
          cache));
  const std::string subfile_path_with_suffix =
      subfile_path_ + section_observation_key.value();

  // Nodegroups
  if constexpr (std::is_same_v<typename ComponentName::chare_type,
                               Parallel::Algorithms::Nodegroup>) {
    auto& local_component = *Parallel::local_branch(
        Parallel::get_parallel_component<ComponentName>(cache));

    // Only want to do this observation once per node for nodegroups. Check if
    // this element is the designated element for this node
    if (designated_node) {
      const double size_in_bytes =
          static_cast<double>(size_of_object_in_bytes(local_component));
      const double size_in_MB = size_in_bytes / 1.0e6;

      Parallel::simple_action<observers::Actions::ContributeReductionData>(
          local_observer,
          observers::ObservationId(observation_value,
                                   subfile_path_with_suffix + ".dat"),
          observers::ArrayComponentId{
              std::add_pointer_t<ParallelComponent>{nullptr},
              Parallel::ArrayIndex<ArrayIndex>(array_index)},
          subfile_path_with_suffix, legend,
          ReductionData{static_cast<double>(observation_value), sys::my_node(),
                        1_st, size_in_MB});
    }
    // Groups
  } else if constexpr (std::is_same_v<typename ComponentName::chare_type,
                                      Parallel::Algorithms::Group>) {
    auto& local_component = *Parallel::local_branch(
        Parallel::get_parallel_component<ComponentName>(cache));

    // Only want to do this observation once per proc for groups. Check if this
    // element is the designated element for this node
    if (designated_proc) {
      const double size_in_bytes =
          static_cast<double>(size_of_object_in_bytes(local_component));
      const double size_in_MB = size_in_bytes / 1.0e6;

      Parallel::simple_action<observers::Actions::ContributeReductionData>(
          local_observer,
          observers::ObservationId(observation_value,
                                   subfile_path_with_suffix + ".dat"),
          observers::ArrayComponentId{
              std::add_pointer_t<ParallelComponent>{nullptr},
              Parallel::ArrayIndex<ArrayIndex>(array_index)},
          subfile_path_with_suffix, legend,
          ReductionData{static_cast<double>(observation_value), sys::my_proc(),
                        1_st, size_in_MB});
    }
    // MutableGlobalCache
  } else if constexpr (std::is_same_v<ComponentName, void>) {
    auto& local_component =
        *Parallel::local_branch(cache.mutable_cache_proxy());
    // Only want to do this observation once per node for the GlobalCache. Check
    // if this element is the designated element for this node
    if (designated_proc) {
      const double size_in_bytes =
          static_cast<double>(size_of_object_in_bytes(local_component));
      const double size_in_MB = size_in_bytes / 1.0e6;

      Parallel::simple_action<observers::Actions::ContributeReductionData>(
          local_observer,
          observers::ObservationId(observation_value,
                                   subfile_path_with_suffix + ".dat"),
          observers::ArrayComponentId{
              std::add_pointer_t<ParallelComponent>{nullptr},
              Parallel::ArrayIndex<ArrayIndex>(array_index)},
          subfile_path_with_suffix, legend,
          ReductionData{static_cast<double>(observation_value), sys::my_proc(),
                        1_st, size_in_MB});
    }
    // Singletons
  } else if constexpr (std::is_same_v<typename ComponentName::chare_type,
                                      Parallel::Algorithms::Singleton>) {
    // Only want one node for this as there is only one singleton so choose node
    // 0
    if (designated_node and sys::my_node() == 0) {
      auto& local_component = *Parallel::local(
          Parallel::get_parallel_component<ComponentName>(cache));
      const double size_in_bytes =
          static_cast<double>(size_of_object_in_bytes(local_component));
      const double size_in_MB = size_in_bytes / 1.0e6;

      auto& observer_writer_proxy = Parallel::get_parallel_component<
          observers::ObserverWriter<Metavariables>>(cache);

      singleton_memory_monitor(
          observer_writer_proxy, subfile_path_with_suffix, legend,
          std::make_tuple(static_cast<double>(observation_value), 0, 1_st,
                          size_in_MB));
    }
  } else if constexpr (std::is_same_v<typename ComponentName::chare_type,
                                      Parallel::Algorithms::Array>) {
    // Now do the elements as all elements will need to contribute to this. We
    // group these by node for simplicity, but all elements will contribute
    auto& local_element =
        *Parallel::local(Parallel::get_parallel_component<ParallelComponent>(
            cache)[array_index]);
    const double size_in_bytes =
        static_cast<double>(size_of_object_in_bytes(local_element));
    const double size_in_MB = size_in_bytes / 1.0e6;

    Parallel::simple_action<observers::Actions::ContributeReductionData>(
        local_observer,
        observers::ObservationId(observation_value,
                                 subfile_path_with_suffix + ".dat"),
        observers::ArrayComponentId{
            std::add_pointer_t<ParallelComponent>{nullptr},
            Parallel::ArrayIndex<ArrayIndex>(array_index)},
        subfile_path_with_suffix, legend,
        ReductionData{static_cast<double>(observation_value), sys::my_node(),
                      1_st, size_in_MB});
  }
}

template <typename ObservationValueTag, typename Metavariables, size_t Index>
void MemoryMonitor<ObservationValueTag, Metavariables, Index>::pup(PUP::er& p) {
  Event::pup(p);
  p | subfile_path_;
}

template <typename ObservationValueTag, typename Metavariables, size_t Index>
PUP::able::PUP_ID
    MemoryMonitor<ObservationValueTag, Metavariables, Index>::my_PUP_ID =
        0;  // NOLINT
/// \endcond

// template <typename ObserveTag, typename Metavariables>
// using memory_monitor_events = tmpl::transform<
//     tmpl::range<size_t, 0,
//                 tmpl::size<typename Metavariables::component_list>::value +
//                 1>,
//     tmpl::bind<MemoryMonitor, tmpl::pin<ObserveTag>,
//     tmpl::pin<Metavariables>,
//                tmpl::_1>>;
}  // namespace Events
