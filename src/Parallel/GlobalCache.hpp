// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines class template GlobalCache.

#pragma once

#include <charm++.h>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <tuple>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataBox/TagTraits.hpp"
#include "Parallel/Callback.hpp"
#include "Parallel/CharmRegistration.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/PupStlCpp17.hpp"
#include "Parallel/Serialize.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/CreateGetTypeAliasOrDefault.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

#include "Parallel/GlobalCache.decl.h"
#include "Parallel/Main.decl.h"

#include "Parallel/Printf.hpp"

/// \cond
template <typename Metavariables>
struct MemoryMonitor;
namespace Parallel {
template <typename Metavariables>
struct GlobalCache;
}  // namespace Parallel
namespace mem_monitor {
template <typename ContributingComponent>
struct ContributeMemoryData;
}  // namespace mem_monitor
/// \endcond

namespace Parallel {

namespace GlobalCache_detail {

template <class GlobalCacheTag, class Metavariables>
using get_matching_tag = typename matching_tag_helper<
    GlobalCacheTag,
    tmpl::append<get_const_global_cache_tags<Metavariables>,
                 get_mutable_global_cache_tags<Metavariables>>>::type;

template <class GlobalCacheTag, class Metavariables>
using type_for_get = typename type_for_get_helper<
    typename get_matching_tag<GlobalCacheTag, Metavariables>::type>::type;

CREATE_GET_TYPE_ALIAS_OR_DEFAULT(component_being_mocked)

template <typename... Tags>
auto make_mutable_cache_tag_storage(tuples::TaggedTuple<Tags...>&& input) {
  return tuples::TaggedTuple<MutableCacheTag<Tags>...>(
      std::make_tuple(std::move(tuples::get<Tags>(input)),
                      std::vector<std::unique_ptr<Callback>>{})...);
}

template <typename ParallelComponent, typename ComponentList>
auto get_component_if_mocked_impl() {
  if constexpr (tmpl::list_contains_v<ComponentList, ParallelComponent>) {
    return ParallelComponent{};
  } else {
    using mock_find = tmpl::find<
        ComponentList,
        std::is_same<get_component_being_mocked_or_default<tmpl::_1, void>,
                     tmpl::pin<ParallelComponent>>>;
    static_assert(
        tmpl::size<mock_find>::value > 0,
        "The requested parallel component (first template argument) is not a "
        "known parallel component (second template argument) or a component "
        "being mocked by one of those components.");
    return tmpl::front<mock_find>{};
  }
}

/// In order to be able to use a mock action testing framework we need to be
/// able to get the correct parallel component from the global cache even when
/// the correct component is a mock. We do this by having the mocked
/// components have a member type alias `component_being_mocked`, and having
/// `Parallel::get_component` check if the component to be retrieved is in the
/// `metavariables::component_list`. If it is not in the `component_list` then
/// we search for a mock component that is mocking the component we are trying
/// to retrieve.
template <typename ComponentList, typename ParallelComponent>
using get_component_if_mocked =
    decltype(get_component_if_mocked_impl<ParallelComponent, ComponentList>());
}  // namespace GlobalCache_detail

/// \cond
template <typename ParallelComponentTag, typename Metavariables>
auto get_parallel_component(GlobalCache<Metavariables>& cache)
    -> Parallel::proxy_from_parallel_component<
        GlobalCache_detail::get_component_if_mocked<
            typename Metavariables::component_list, ParallelComponentTag>>&;

template <typename ParallelComponentTag, typename Metavariables>
auto get_parallel_component(const GlobalCache<Metavariables>& cache)
    -> const Parallel::proxy_from_parallel_component<
        GlobalCache_detail::get_component_if_mocked<
            typename Metavariables::component_list, ParallelComponentTag>>&;
/// \endcond

/// \ingroup ParallelGroup
/// A Charm++ chare that caches mutable data once per Charm++ core.
///
/// `MutableGlobalCache` is not intended to be visible to the end user; its
/// interface is via the `GlobalCache` member functions
/// `mutable_cache_item_is_ready`, `mutate`, and `get`.
/// Accordingly, most documentation of `MutableGlobalCache` is provided
/// in the relevant `GlobalCache` member functions.
template <typename Metavariables>
class MutableGlobalCache : public CBase_MutableGlobalCache<Metavariables> {
 public:
  // Even though the MutableGlobalCache doesn't run the Algorithm, this type
  // alias helps in identifying that the MutableGlobalCache is a Group using
  // Parallel::is_group_v
  using chare_type = Parallel::Algorithms::Group;

  explicit MutableGlobalCache(tuples::tagged_tuple_from_typelist<
                              get_mutable_global_cache_tags<Metavariables>>
                                  mutable_global_cache);
  explicit MutableGlobalCache(CkMigrateMessage* msg)
      : CBase_MutableGlobalCache<Metavariables>(msg) {}

  ~MutableGlobalCache() override {
    (void)Parallel::charmxx::RegisterChare<
        MutableGlobalCache<Metavariables>,
        CkIndex_MutableGlobalCache<Metavariables>>::registrar;
  }
  /// \cond
  MutableGlobalCache() = default;
  MutableGlobalCache(const MutableGlobalCache&) = default;
  MutableGlobalCache& operator=(const MutableGlobalCache&) = default;
  MutableGlobalCache(MutableGlobalCache&&) = default;
  MutableGlobalCache& operator=(MutableGlobalCache&&) = default;
  /// \endcond

  template <typename GlobalCacheTag>
  auto get() const
      -> const GlobalCache_detail::type_for_get<GlobalCacheTag, Metavariables>&;

  // Entry method to mutate the object identified by `GlobalCacheTag`.
  // Internally calls Function::apply(), where
  // Function is a struct, and Function::apply is a user-defined
  // static function that mutates the object.  Function::apply() takes
  // as its first argument a gsl::not_null pointer to the object named
  // by the GlobalCacheTag, and then the contents of 'args' as
  // subsequent arguments.  Called via `GlobalCache::mutate`.
  template <typename GlobalCacheTag, typename Function, typename... Args>
  void mutate(const std::tuple<Args...>& args);

  // Entry method that computes the size of the local branch of the
  // MutableGlobalCache and sends it to the MemoryMonitor parallel component.
  //
  // Note: This can only be called if the MemoryMonitor component is in the
  // `component_list` of the metavariables. Also can't be called in the testing
  // framework. Trying to do either of these will result in an ERROR.
  void compute_size_for_memory_monitor(
      CProxy_GlobalCache<Metavariables> global_cache_proxy, const double time);

  // Not an entry method, and intended to be called only from
  // `GlobalCache` via the free function
  // `Parallel::mutable_cache_item_is_ready`.  See the free function
  // `Parallel::mutable_cache_item_is_ready` for documentation.
  template <typename GlobalCacheTag, typename Function>
  bool mutable_cache_item_is_ready(const Function& function);

  void pup(PUP::er& p) override;  // NOLINT

 private:
  tuples::tagged_tuple_from_typelist<
      get_mutable_global_cache_tag_storage<Metavariables>>
      mutable_global_cache_{};
};

template <typename Metavariables>
MutableGlobalCache<Metavariables>::MutableGlobalCache(
    tuples::tagged_tuple_from_typelist<
        get_mutable_global_cache_tags<Metavariables>>
        mutable_global_cache)
    : mutable_global_cache_(GlobalCache_detail::make_mutable_cache_tag_storage(
          std::move(mutable_global_cache))) {}

template <typename Metavariables>
template <typename GlobalCacheTag>
auto MutableGlobalCache<Metavariables>::get() const
    -> const GlobalCache_detail::type_for_get<GlobalCacheTag, Metavariables>& {
  using tag = MutableCacheTag<
      GlobalCache_detail::get_matching_tag<GlobalCacheTag, Metavariables>>;
  if constexpr (tt::is_a_v<std::unique_ptr, typename tag::tag::type>) {
    return *(std::get<0>(tuples::get<tag>(mutable_global_cache_)));
  } else {
    return std::get<0>(tuples::get<tag>(mutable_global_cache_));
  }
}

template <typename Metavariables>
template <typename GlobalCacheTag, typename Function>
bool MutableGlobalCache<Metavariables>::mutable_cache_item_is_ready(
    const Function& function) {
  using tag = MutableCacheTag<GlobalCache_detail::get_matching_mutable_tag<
      GlobalCacheTag, Metavariables>>;
  std::unique_ptr<Callback> optional_callback{};
  if constexpr (tt::is_a_v<std::unique_ptr, typename tag::tag::type>) {
    optional_callback =
        function(*(std::get<0>(tuples::get<tag>(mutable_global_cache_))));
  } else {
    optional_callback =
        function(std::get<0>(tuples::get<tag>(mutable_global_cache_)));
  }
  if (optional_callback) {
    std::get<1>(tuples::get<tag>(mutable_global_cache_))
        .push_back(std::move(optional_callback));
    if (std::get<1>(tuples::get<tag>(mutable_global_cache_)).size() > 20000) {
      ERROR("The number of callbacks in MutableGlobalCache for tag "
            << pretty_type::short_name<GlobalCacheTag>()
            << " has gotten too large, and may be growing without bound");
    }
    return false;
  } else {
    // The user-defined `function` didn't specify a callback, which
    // means that the item is ready.
    return true;
  }
}

template <typename Metavariables>
template <typename GlobalCacheTag, typename Function, typename... Args>
void MutableGlobalCache<Metavariables>::mutate(
    const std::tuple<Args...>& args) {
  (void)Parallel::charmxx::RegisterMutableGlobalCacheMutate<
      Metavariables, GlobalCacheTag, Function, Args...>::registrar;
  using tag = MutableCacheTag<GlobalCache_detail::get_matching_mutable_tag<
      GlobalCacheTag, Metavariables>>;

  // Do the mutate.
  std::apply(
      [this](const auto&... local_args) {
        Function::apply(make_not_null(&std::get<0>(
                            tuples::get<tag>(mutable_global_cache_))),
                        local_args...);
      },
      args);

  // A callback might call mutable_cache_item_is_ready, which might
  // add yet another callback to the vector of callbacks.  We don't
  // want to immediately invoke this new callback and we don't want to
  // remove it from the vector of callbacks before it is
  // invoked. Therefore, std::move the vector of callbacks into a
  // temporary vector, clear the original vector, and invoke the callbacks
  // in the temporary vector.
  std::vector<std::unique_ptr<Callback>> callbacks(
      std::move(std::get<1>(tuples::get<tag>(mutable_global_cache_))));
  std::get<1>(tuples::get<tag>(mutable_global_cache_)).clear();
  std::get<1>(tuples::get<tag>(mutable_global_cache_)).shrink_to_fit();

  // Invoke the callbacks.  Any new callbacks that are added to the
  // list (if a callback calls mutable_cache_item_is_ready) will be
  // saved and will not be invoked here.
  for (auto& callback : callbacks) {
    callback->invoke();
  }
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif  // defined(__GNUC__) && !defined(__clang__)
template <typename Metavariables>
void MutableGlobalCache<Metavariables>::compute_size_for_memory_monitor(
    CProxy_GlobalCache<Metavariables> global_cache_proxy, const double time) {
  if constexpr (tmpl::list_contains_v<typename Metavariables::component_list,
                                      MemoryMonitor<Metavariables>>) {
    const double size_in_bytes =
        static_cast<double>(size_of_object_in_bytes(*this));
    const double size_in_MB = size_in_bytes / 1.0e6;

    auto& cache = *Parallel::local_branch(global_cache_proxy);

    auto& mem_monitor_proxy =
        Parallel::get_parallel_component<MemoryMonitor<Metavariables>>(cache);

    // We use sys::my_proc() here rather than
    // Parallel::my_proc(*Parallel::local_branch(this->thisProxy)) because the
    // MutableGlobalCache doesn't have a `my_proc()` member to call
    const int my_proc = sys::my_proc();

    Parallel::simple_action<
        mem_monitor::ContributeMemoryData<MutableGlobalCache<Metavariables>>>(
        mem_monitor_proxy, time, my_proc, size_in_MB);
  } else {
    (void)global_cache_proxy;
    (void)time;
    ERROR(
        "This entry method of the MutableGlobalCache can only be called if the "
        "MemoryMonitor is in the component list in the metavariables.\n");
  }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__)

template <typename Metavariables>
void MutableGlobalCache<Metavariables>::pup(PUP::er& p) {
  p | mutable_global_cache_;
}

/// \ingroup ParallelGroup
/// A Charm++ chare that caches constant data once per Charm++ node or
/// non-constant data once per Charm++ core.
///
/// `Metavariables` must define the following metavariables:
///   - `component_list`   typelist of ParallelComponents
///   - `const_global_cache_tags`   (possibly empty) typelist of tags of
///     constant data
///   - `mutable_global_cache_tags` (possibly empty) typelist of tags of
///     non-constant data
///
/// The tag lists for the const items added to the GlobalCache is created by
/// combining the following tag lists:
///   - `Metavariables::const_global_cache_tags` which should contain only those
///     tags that cannot be added from the other tag lists below.
///   - `Component::const_global_cache_tags` for each `Component` in
///     `Metavariables::component_list` which should contain the tags needed by
///     any simple actions called on the Component, as well as tags need by the
///     `allocate_array` function of an array component.  The type alias may be
///     omitted for an empty list.
///   - `Action::const_global_cache_tags` for each `Action` in the
///     `phase_dependent_action_list` of each `Component` of
///     `Metavariables::component_list` which should contain the tags needed by
///     that  Action.  The type alias may be omitted for an empty list.
///
/// The tag lists for the non-const items added to the GlobalCache is created
/// by combining exactly the same tag lists as for the const items, except with
/// `const_global_cache_tags` replaced by `mutable_global_cache_tags`.
///
/// The tags in the `const_global_cache_tags` and
/// `mutable_global_cache_tags` type lists are db::SimpleTag%s that
/// have a `using option_tags` type alias and a static function
/// `create_from_options` that are used to create the constant data
/// from input file options.
///
/// References to const items in the GlobalCache are also added to the
/// db::DataBox of each `Component` in the
/// `Metavariables::component_list` with the same tag with which they
/// were inserted into the GlobalCache.  References to non-const items
/// in the GlobalCache are not added to the db::DataBox.
template <typename Metavariables>
class GlobalCache : public CBase_GlobalCache<Metavariables> {
  using parallel_component_tag_list = tmpl::transform<
      typename Metavariables::component_list,
      tmpl::bind<
          tmpl::type_,
          tmpl::bind<Parallel::proxy_from_parallel_component, tmpl::_1>>>;

 public:
  using proxy_type = CProxy_GlobalCache<Metavariables>;
  using main_proxy_type = CProxy_Main<Metavariables>;
  using mutable_cache_proxy_type = CProxy_MutableGlobalCache<Metavariables>;
  /// Access to the Metavariables template parameter
  using metavariables = Metavariables;
  /// Typelist of the ParallelComponents stored in the GlobalCache
  using component_list = typename Metavariables::component_list;
  // Even though the GlobalCache doesn't run the Algorithm, this type alias
  // helps in identifying that the GlobalCache is a Nodegroup using
  // Parallel::is_nodegroup_v
  using chare_type = Parallel::Algorithms::Nodegroup;

  /// Constructor used only by the ActionTesting framework and other
  /// non-charm++ tests that don't know about proxies.
  GlobalCache(tuples::tagged_tuple_from_typelist<
                  get_const_global_cache_tags<Metavariables>>
                  const_global_cache,
              MutableGlobalCache<Metavariables>* mutable_global_cache,
              std::optional<main_proxy_type> main_proxy = std::nullopt);

  /// Constructor used by Main and anything else that is charm++ aware.
  GlobalCache(
      tuples::tagged_tuple_from_typelist<
          get_const_global_cache_tags<Metavariables>>
          const_global_cache,
      CProxy_MutableGlobalCache<Metavariables> mutable_global_cache_proxy,
      std::optional<main_proxy_type> main_proxy = std::nullopt);

  explicit GlobalCache(CkMigrateMessage* msg)
      : CBase_GlobalCache<Metavariables>(msg) {}

  ~GlobalCache() override {
    (void)Parallel::charmxx::RegisterChare<
        GlobalCache<Metavariables>,
        CkIndex_GlobalCache<Metavariables>>::registrar;
  }
  /// \cond
  GlobalCache() = default;
  GlobalCache(const GlobalCache&) = default;
  GlobalCache& operator=(const GlobalCache&) = default;
  GlobalCache(GlobalCache&&) = default;
  GlobalCache& operator=(GlobalCache&&) = default;
  /// \endcond

  /// Entry method to set the ParallelComponents (should only be called once)
  void set_parallel_components(
      tuples::tagged_tuple_from_typelist<parallel_component_tag_list>&&
          parallel_components,
      const CkCallback& callback);

  /// Returns whether the object referred to by `GlobalCacheTag`
  /// (which must be a mutable cache tag) is ready to be accessed by a
  /// `get` call.
  ///
  /// `function` is a user-defined invokable that:
  /// - takes one argument: a const reference to the object referred to by the
  ///   `GlobalCacheTag`.
  /// - if the data is ready, returns a default constructed
  ///   `std::unique_ptr<CallBack>`
  /// - if the data is not ready, returns a `std::unique_ptr<CallBack>`,
  ///   where the `Callback` will re-invoke the current action on the
  ///   current parallel component.
  template <typename GlobalCacheTag, typename Function>
  bool mutable_cache_item_is_ready(const Function& function);

  /// Mutates the non-const object identified by GlobalCacheTag.
  /// \requires `GlobalCacheTag` is a tag in `mutable_global_cache_tags`
  /// defined by the Metavariables and in Actions.
  ///
  /// Internally calls `Function::apply()`, where `Function` is a
  /// user-defined struct and `Function::apply()` is a user-defined
  /// static function that mutates the object.  `Function::apply()`
  /// takes as its first argument a `gsl::not_null` pointer to the
  /// object named by the GlobalCacheTag, and takes the contents of
  /// `args` as subsequent arguments.
  template <typename GlobalCacheTag, typename Function, typename... Args>
  void mutate(const std::tuple<Args...>& args);

  /// Entry method that computes the size of the local branch of the
  /// GlobalCache and sends it to the MemoryMonitor parallel component.
  ///
  /// \note This can only be called if the MemoryMonitor component is in the
  /// `component_list` of the metavariables. Also can't be called in the testing
  /// framework. Trying to do either of these will result in an ERROR.
  void compute_size_for_memory_monitor(const double time);

  /// Retrieve the proxy to the global cache
  proxy_type get_this_proxy();

  void pup(PUP::er& p) override;  // NOLINT

  /// Retrieve the proxy to the Main chare (or std::nullopt if the proxy has not
  /// been set).
  std::optional<main_proxy_type> get_main_proxy();

  /// Returns `true` if the GlobalCache was constructed with a Proxy to the
  /// MutableGlobalCache (i.e. the GlobalCache is charm-aware).
  bool mutable_global_cache_proxy_is_set() const;

  /// Returns a proxy to the MutableGlobalCache if the proxy has been set. If it
  /// hasn't been set, an ERROR will occur meaning that this method can't be
  /// used in the testing framework.
  mutable_cache_proxy_type mutable_cache_proxy();

 private:
  // clang-tidy: false positive, redundant declaration
  template <typename GlobalCacheTag, typename MV>
  friend auto get(const GlobalCache<MV>& cache)  // NOLINT
      -> const GlobalCache_detail::type_for_get<GlobalCacheTag, MV>&;

  // clang-tidy: false positive, redundant declaration
  template <typename ParallelComponentTag, typename MV>
  friend auto get_parallel_component(  // NOLINT
      GlobalCache<MV>& cache)
      -> Parallel::proxy_from_parallel_component<
          GlobalCache_detail::get_component_if_mocked<
              typename MV::component_list, ParallelComponentTag>>&;

  // clang-tidy: false positive, redundant declaration
  template <typename ParallelComponentTag, typename MV>
  friend auto get_parallel_component(  // NOLINT
      const GlobalCache<MV>& cache)
      -> const Parallel::proxy_from_parallel_component<
          GlobalCache_detail::get_component_if_mocked<
              typename MV::component_list,
              ParallelComponentTag>>&;  // NOLINT

  tuples::tagged_tuple_from_typelist<get_const_global_cache_tags<Metavariables>>
      const_global_cache_{};
  tuples::tagged_tuple_from_typelist<parallel_component_tag_list>
      parallel_components_{};
  // We store both a pointer and a proxy to the MutableGlobalCache.
  // There is both a pointer and a proxy because we want to use
  // MutableGlobalCache in production (where it must be charm-aware)
  // and for simple testing (which we want to do in a non-charm-aware
  // context for simplicity).  If the charm-aware constructor is
  // used, then the pointer is set to nullptr and the proxy is set.
  // If the non-charm-aware constructor is used, the the pointer is
  // set and the proxy is ignored.  The member functions that need the
  // MutableGlobalCache should use the pointer if it is not nullptr,
  // otherwise use the proxy.
  MutableGlobalCache<Metavariables>* mutable_global_cache_{nullptr};
  mutable_cache_proxy_type mutable_global_cache_proxy_{};
  bool parallel_components_have_been_set_{false};
  std::optional<main_proxy_type> main_proxy_;
};

template <typename Metavariables>
GlobalCache<Metavariables>::GlobalCache(
    tuples::tagged_tuple_from_typelist<
        get_const_global_cache_tags<Metavariables>>
        const_global_cache,
    MutableGlobalCache<Metavariables>* mutable_global_cache,
    std::optional<main_proxy_type> main_proxy)
    : const_global_cache_(std::move(const_global_cache)),
      mutable_global_cache_(mutable_global_cache),
      main_proxy_(std::move(main_proxy)) {
  ASSERT(mutable_global_cache_ != nullptr,
         "GlobalCache: Do not construct with a nullptr!");
}

template <typename Metavariables>
GlobalCache<Metavariables>::GlobalCache(
    tuples::tagged_tuple_from_typelist<
        get_const_global_cache_tags<Metavariables>>
        const_global_cache,
    CProxy_MutableGlobalCache<Metavariables> mutable_global_cache_proxy,
    std::optional<main_proxy_type> main_proxy)
    : const_global_cache_(std::move(const_global_cache)),
      mutable_global_cache_(nullptr),
      mutable_global_cache_proxy_(std::move(mutable_global_cache_proxy)),
      main_proxy_(std::move(main_proxy)) {}

template <typename Metavariables>
void GlobalCache<Metavariables>::set_parallel_components(
    tuples::tagged_tuple_from_typelist<parallel_component_tag_list>&&
        parallel_components,
    const CkCallback& callback) {
  ASSERT(!parallel_components_have_been_set_,
         "Can only set the parallel_components once");
  parallel_components_ = std::move(parallel_components);
  parallel_components_have_been_set_ = true;
  this->contribute(callback);
}

template <typename Metavariables>
template <typename GlobalCacheTag, typename Function>
bool GlobalCache<Metavariables>::mutable_cache_item_is_ready(
    const Function& function) {
  if (mutable_global_cache_proxy_is_set()) {
    return Parallel::local_branch(mutable_global_cache_proxy_)
        ->template mutable_cache_item_is_ready<GlobalCacheTag>(function);
  } else {
    return mutable_global_cache_
        ->template mutable_cache_item_is_ready<GlobalCacheTag>(function);
  }
}

template <typename Metavariables>
template <typename GlobalCacheTag, typename Function, typename... Args>
void GlobalCache<Metavariables>::mutate(const std::tuple<Args...>& args) {
  (void)Parallel::charmxx::RegisterGlobalCacheMutate<
      Metavariables, GlobalCacheTag, Function, Args...>::registrar;
  if (mutable_global_cache_proxy_is_set()) {
    // charm-aware version: Mutate the variable on all PEs on this node.
    for (auto pe = CkNodeFirst(CkMyNode());
         pe < CkNodeFirst(CkMyNode()) + CkNodeSize(CkMyNode()); ++pe) {
      mutable_global_cache_proxy_[pe].template mutate<GlobalCacheTag, Function>(
          args);
    }
  } else {
    // version that bypasses proxies.  Just call the function.
    mutable_global_cache_->template mutate<GlobalCacheTag, Function>(args);
  }
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif  // defined(__GNUC__) && !defined(__clang__)
template <typename Metavariables>
void GlobalCache<Metavariables>::compute_size_for_memory_monitor(
    const double time) {
  if constexpr (tmpl::list_contains_v<typename Metavariables::component_list,
                                      MemoryMonitor<Metavariables>>) {
    const double size_in_bytes =
        static_cast<double>(size_of_object_in_bytes(*this));
    const double size_in_MB = size_in_bytes / 1.0e6;

    auto& mem_monitor_proxy =
        Parallel::get_parallel_component<MemoryMonitor<Metavariables>>(*this);

    // We use sys::my_node() here rather than
    // Parallel::my_node(*Parallel::local_branch(this->get_this_proxy()))
    // because the MutableGlobalCache doesn't have a `my_node()` member to call
    const int my_node = sys::my_node();

    Parallel::simple_action<
        mem_monitor::ContributeMemoryData<GlobalCache<Metavariables>>>(
        mem_monitor_proxy, time, my_node, size_in_MB);
  } else {
    (void)time;
    ERROR(
        "This entry method of the GlobalCache can only be called if the "
        "MemoryMonitor is in the component list in the metavariables.\n");
  }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__)

template <typename Metavariables>
typename Parallel::GlobalCache<Metavariables>::proxy_type
GlobalCache<Metavariables>::get_this_proxy() {
  return this->thisProxy;
}

template <typename Metavariables>
std::optional<typename Parallel::GlobalCache<Metavariables>::main_proxy_type>
GlobalCache<Metavariables>::get_main_proxy() {
  if (main_proxy_.has_value()) {
    return main_proxy_;
  } else {
    ERROR(
        "Attempting to retrieve the main proxy in a context in which the main "
        "proxy has not been supplied to the constructor.");
  }
}

template <typename Metavariables>
bool GlobalCache<Metavariables>::mutable_global_cache_proxy_is_set() const {
  return mutable_global_cache_ == nullptr;
}

template <typename Metavariables>
typename Parallel::GlobalCache<Metavariables>::mutable_cache_proxy_type
GlobalCache<Metavariables>::mutable_cache_proxy() {
  if (not mutable_global_cache_proxy_is_set()) {
    ERROR(
        "Cannot get a proxy to the mutable global cache because the proxy "
        "isn't set.");
  }
  return mutable_global_cache_proxy_;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif  // defined(__GNUC__) && !defined(__clang__)
template <typename Metavariables>
void GlobalCache<Metavariables>::pup(PUP::er& p) {
  p | const_global_cache_;
  p | parallel_components_;
  p | mutable_global_cache_proxy_;
  p | main_proxy_;
  p | parallel_components_have_been_set_;
  if (not p.isUnpacking() and mutable_global_cache_ != nullptr) {
    ERROR(
        "Cannot serialize the const global cache when the mutable global cache "
        "is set to a local pointer. If this occurs in a unit test, avoid the "
        "serialization. If this occurs in a production executable, be sure "
        "that the MutableGlobalCache is accessed by a charm proxy.");
  }
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__)

/// @{
/// \ingroup ParallelGroup
/// \brief Access the Charm++ proxy associated with a ParallelComponent
///
/// \requires ParallelComponentTag is a tag in component_list
///
/// \returns a Charm++ proxy that can be used to call an entry method on the
/// chare(s)
template <typename ParallelComponentTag, typename Metavariables>
auto get_parallel_component(GlobalCache<Metavariables>& cache)
    -> Parallel::proxy_from_parallel_component<
        GlobalCache_detail::get_component_if_mocked<
            typename Metavariables::component_list, ParallelComponentTag>>& {
  return tuples::get<tmpl::type_<Parallel::proxy_from_parallel_component<
      GlobalCache_detail::get_component_if_mocked<
          typename Metavariables::component_list, ParallelComponentTag>>>>(
      cache.parallel_components_);
}

template <typename ParallelComponentTag, typename Metavariables>
auto get_parallel_component(const GlobalCache<Metavariables>& cache)
    -> const Parallel::proxy_from_parallel_component<
        GlobalCache_detail::get_component_if_mocked<
            typename Metavariables::component_list, ParallelComponentTag>>& {
  return tuples::get<tmpl::type_<Parallel::proxy_from_parallel_component<
      GlobalCache_detail::get_component_if_mocked<
          typename Metavariables::component_list, ParallelComponentTag>>>>(
      cache.parallel_components_);
}
/// @}

/// @{
/// \ingroup ParallelGroup
/// \brief Access data in the cache
///
/// \requires GlobalCacheTag is a tag in the `mutable_global_cache_tags`
/// or `const_global_cache_tags` defined by the Metavariables and in Actions.
///
/// \returns a constant reference to an object in the cache
template <typename GlobalCacheTag, typename Metavariables>
auto get(const GlobalCache<Metavariables>& cache)
    -> const GlobalCache_detail::type_for_get<GlobalCacheTag, Metavariables>& {
  // We check if the tag is to be retrieved directly or via a base class
  using tag =
      GlobalCache_detail::get_matching_tag<GlobalCacheTag, Metavariables>;
  if constexpr (is_in_mutable_global_cache<Metavariables, GlobalCacheTag>) {
    // Tag is not in the const tags, so use MutableGlobalCache
    if (cache.mutable_global_cache_proxy_is_set()) {
      const auto& local_mutable_cache =
          *Parallel::local_branch(cache.mutable_global_cache_proxy_);
      return local_mutable_cache.template get<GlobalCacheTag>();
    } else {
      return cache.mutable_global_cache_->template get<GlobalCacheTag>();
    }
  } else {
    // Tag is in the const tags, so use const_global_cache_
    if constexpr (tt::is_a_v<std::unique_ptr, typename tag::type>) {
      return *(tuples::get<tag>(cache.const_global_cache_));
    } else {
      return tuples::get<tag>(cache.const_global_cache_);
    }
  }
}

/// \ingroup ParallelGroup
/// \brief Returns whether the object identified by `GlobalCacheTag`
/// is ready to be accessed by `get`.
///
/// \requires `GlobalCacheTag` is a tag in `mutable_global_cache_tags`
/// defined by the Metavariables and in Actions.
///
/// \requires `function` is a user-defined invokable that takes one argument:
/// a const reference to the object referred to by the
/// `GlobalCacheTag`.  `function` returns a
/// `std::unique_ptr<CallBack>` that determines the readiness. To
/// indicate that the item is ready, the `std::unique_ptr` returned
/// by `function` must be nullptr; in this case
/// `mutable_cache_item_is_ready` returns true. To indicate that the
/// item is not ready, the `std::unique_ptr` returned by `function`
/// must be valid; in this case, `mutable_cache_item_is_ready`
/// appends the `std::unique_ptr<Callback>` to an
/// internal list of callbacks to be called on `mutate`, and then
/// returns false.
template <typename GlobalCacheTag, typename Function, typename Metavariables>
bool mutable_cache_item_is_ready(GlobalCache<Metavariables>& cache,
                                 const Function& function) {
  return cache.template mutable_cache_item_is_ready<GlobalCacheTag>(function);
}

/// \ingroup ParallelGroup
///
/// \brief Mutates non-const data in the cache, by calling `Function::apply()`
///
/// \requires `GlobalCacheTag` is a tag in tag_list.
/// \requires `Function` is a struct with a static void `apply()`
/// function that mutates the object. `Function::apply()` takes as its
/// first argument a `gsl::not_null` pointer to the object named by
/// the `GlobalCacheTag`, and takes `args` as
/// subsequent arguments.
template <typename GlobalCacheTag, typename Function, typename Metavariables,
          typename... Args>
void mutate(GlobalCache<Metavariables>& cache, Args&&... args) {
  if (cache.mutable_global_cache_proxy_is_set()) {
    cache.thisProxy.template mutate<GlobalCacheTag, Function>(
        std::make_tuple<Args...>(std::forward<Args>(args)...));
  } else {
    cache.template mutate<GlobalCacheTag, Function>(
        std::make_tuple<Args...>(std::forward<Args>(args)...));
  }
}

namespace Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup ParallelGroup
/// Tag to retrieve the `Parallel::GlobalCache` from the DataBox.
struct GlobalCache : db::BaseTag {};

template <class Metavariables>
struct GlobalCacheProxy : db::SimpleTag {
  using type = CProxy_GlobalCache<Metavariables>;
};

template <class Metavariables>
struct GlobalCacheImpl : GlobalCache, db::SimpleTag {
  using type = Parallel::GlobalCache<Metavariables>*;
  static std::string name() { return "GlobalCache"; }
};

template <class Metavariables>
struct GlobalCacheImplCompute : GlobalCacheImpl<Metavariables>, db::ComputeTag {
  using base = GlobalCacheImpl<Metavariables>;
  using argument_tags = tmpl::list<GlobalCacheProxy<Metavariables>>;
  using return_type = Parallel::GlobalCache<Metavariables>*;
  static void function(
      const gsl::not_null<Parallel::GlobalCache<Metavariables>**>
          local_branch_of_global_cache,
      const CProxy_GlobalCache<Metavariables>& global_cache_proxy) {
    *local_branch_of_global_cache = Parallel::local_branch(global_cache_proxy);
  }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ParallelGroup
/// Tag used to retrieve data from the `Parallel::GlobalCache`. This is the
/// recommended way for compute tags to retrieve data out of the global cache.
template <class CacheTag>
struct FromGlobalCache : CacheTag, db::ReferenceTag {
  static_assert(db::is_simple_tag_v<CacheTag>);
  using base = CacheTag;
  using parent_tag = GlobalCache;

  template <class Metavariables>
  static const auto& get(
      const Parallel::GlobalCache<Metavariables>* const& cache) {
    return Parallel::get<CacheTag>(*cache);
  }

  using argument_tags = tmpl::list<parent_tag>;
};
}  // namespace Tags
}  // namespace Parallel

namespace PUP {
/// \cond
// Warning! This is an invalid and kludgey pupper, because when unpacking it
// produces a `nullptr` rather than a valid `GlobalCache*`. This function is
// provided _only_ to enable putting a `GlobalCache*` in the DataBox.
//
// SpECTRE parallel components with a `GlobalCache*` in their DataBox should
// set this pointer using a compute item that calls `Parallel::local_branch` on
// a stored proxy. When deserializing the DataBox, these components should call
// `db:mutate<GlobalCacheProxy>(box)` to force the DataBox to update the
// pointer from the new Charm++ proxy.
//
// We do not currently anticipate needing to (de)serialize a `GlobalCache*`
// outside of a DataBox. But if this need arises, it will be necessary to
// provide a non-kludgey pupper here.
//
// Correctly (de)serializing the `GlobalCache*` would require obtaining a
// `CProxy_GlobalCache` item and calling `Parallel::local_branch` on it ---
// just as in the workflow described above, but within the pupper vs in the
// DataBox. But this strategy fails when restarting from a checkpoint file,
// because calling `Parallel::local_branch` may not be well-defined in the
// unpacking pup call when all Charm++ components may not yet be fully restored.
// This difficulty is why we instead write an invalid pupper here.
//
// In future versions of Charm++, the pup function may know whether it is
// called in the context of checkpointing, load balancing, etc. This knowledge
// would enable us to write a valid pupper for non-checkpoint contexts, and
// return a `nullptr` only when restoring from checkpoint.
template <typename Metavariables>
inline void operator|(PUP::er& p,  // NOLINT
                      Parallel::GlobalCache<Metavariables>*& t) {
  if (p.isUnpacking()) {
    t = nullptr;
  }
}
/// \endcond
}  // namespace PUP

#define CK_TEMPLATES_ONLY
#include "Parallel/GlobalCache.def.h"
#undef CK_TEMPLATES_ONLY
