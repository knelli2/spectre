// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines the Charm++ mainchare.

#pragma once

#include <boost/program_options.hpp>
#include <charm++.h>
#include <initializer_list>
#include <pup.h>
#include <regex>
#include <string>
#include <type_traits>
#include <unordered_set>

#include "Informer/Informer.hpp"
#include "Options/ParseOptions.hpp"
#include "Parallel/AlgorithmMetafunctions.hpp"
#include "Parallel/CharmRegistration.hpp"
#include "Parallel/CreateFromOptions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/PhaseControl/InitializePhaseChangeDecisionData.hpp"
#include "Parallel/PhaseControl/PhaseControlTags.hpp"
#include "Parallel/PhaseControlReductionHelpers.hpp"
#include "Parallel/Printf.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/Tags/ResourceInfo.hpp"
#include "Parallel/TypeTraits.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Formaline.hpp"
#include "Utilities/Overloader.hpp"
#include "Utilities/System/Exit.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/CreateGetTypeAliasOrDefault.hpp"

#include "Parallel/Main.decl.h"

namespace Parallel {
/// \ingroup ParallelGroup
/// The main function of a Charm++ executable.
/// See [the Parallelization documentation](group__ParallelGroup.html#details)
/// for an overview of Metavariables, Phases, and parallel components.
template <typename Metavariables>
class Main : public CBase_Main<Metavariables> {
 public:
  using component_list = typename Metavariables::component_list;
  using const_global_cache_tags = get_const_global_cache_tags<Metavariables>;
  using mutable_global_cache_tags =
      get_mutable_global_cache_tags<Metavariables>;

  using phase_change_tags_and_combines_list =
      PhaseControl::get_phase_change_tags<Metavariables>;
  /// \cond HIDDEN_SYMBOLS
  /// The constructor used to register the class
  explicit Main(const Parallel::charmxx::
                    MainChareRegistrationConstructor& /*used_for_reg*/) {}
  ~Main() override {
    (void)Parallel::charmxx::RegisterChare<
        Main<Metavariables>, CkIndex_Main<Metavariables>>::registrar;
  }
  Main(const Main&) = default;
  Main& operator=(const Main&) = default;
  Main(Main&&) = default;
  Main& operator=(Main&&) = default;
  /// \endcond

  explicit Main(CkArgMsg* msg);
  explicit Main(CkMigrateMessage* msg);

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override;

  /// Allocate the initial elements of array components, and then execute the
  /// initialization phase on each component
  void allocate_array_components_and_execute_initialization_phase();

  /// Determine the next phase of the simulation and execute it.
  void execute_next_phase();

  /// Place the Charm++ call that starts load balancing
  ///
  /// \details This call is wrapped within an entry method so that it may be
  /// used as the callback after a quiescence detection.
  void start_load_balance();

  /// Place the Charm++ call that starts writing a checkpoint
  ///
  /// \details This call is wrapped within an entry method so that it may be
  /// used as the callback after a quiescence detection.
  void start_write_checkpoint();

  /// Reduction target for data used in phase change decisions.
  ///
  /// It is required that the `Parallel::ReductionData` holds a single
  /// `tuples::TaggedTuple`.
  template <typename InvokeCombine, typename... Tags>
  void phase_change_reduction(
      ReductionData<ReductionDatum<tuples::TaggedTuple<Tags...>, InvokeCombine,
                                   funcl::Identity, std::index_sequence<>>>
          reduction_data);

 private:
  // Allocate all singleton components and return procs to ignore
  std::unordered_set<size_t> allocate_singleton_components(bool avoid_proc_0);

  // Return the dir name for the next Charm++ checkpoint as well as the pieces
  // from which the name is built up: the basename and the padding. This is a
  // "detail" function so that these pieces can be defined in one place only.
  std::tuple<std::string, std::string, size_t> checkpoint_dir_basename_pad()
      const;

  // Return the dir name for the next Charm++ checkpoint; check and error if
  // this name already exists and writing the checkpoint would be unsafe.
  std::string checkpoint_dir() const;

  // Check if future checkpoint dirs are available; error if any already exist.
  void check_future_checkpoint_dirs_available() const;

  template <typename ParallelComponent>
  using parallel_component_options =
      Parallel::get_option_tags<typename ParallelComponent::initialization_tags,
                                Metavariables>;
  using option_list = tmpl::remove_duplicates<tmpl::flatten<tmpl::list<
      Parallel::get_option_tags<const_global_cache_tags, Metavariables>,
      Parallel::get_option_tags<mutable_global_cache_tags, Metavariables>,
      tmpl::transform<component_list,
                      tmpl::bind<parallel_component_options, tmpl::_1>>>>>;
  using parallel_component_tag_list = tmpl::transform<
      component_list,
      tmpl::bind<
          tmpl::type_,
          tmpl::bind<Parallel::proxy_from_parallel_component, tmpl::_1>>>;

  typename Metavariables::Phase current_phase_{
      Metavariables::Phase::Initialization};
  CProxy_MutableGlobalCache<Metavariables> mutable_global_cache_proxy_;
  CProxy_GlobalCache<Metavariables> global_cache_proxy_;
  detail::CProxy_AtSyncIndicator<Metavariables> at_sync_indicator_proxy_;
  // This is only used during startup, and will be cleared after all
  // the chares are created.  It is a member variable because passing
  // local state through charm callbacks is painful.
  tuples::tagged_tuple_from_typelist<option_list> options_{};
  // type to be determined by the collection of available phase changers in the
  // Metavariables
  tuples::tagged_tuple_from_typelist<phase_change_tags_and_combines_list>
      phase_change_decision_data_;
  size_t checkpoint_dir_counter_ = 0_st;
};

namespace detail {

// Charm++ AtSync effectively requires an additional global sync to the
// quiescence detection we do for switching phases. However, AtSync only needs
// to be called for one array to trigger the sync-based load balancing, so the
// AtSyncIndicator is a silly hack to have a centralized indication to start
// load-balancing. It participates in the `AtSync` barrier, but is not
// migratable, and should only be constructed by the `Main` chare on the same
// processor as the `Main` chare. When load-balancing occurs, main invokes the
// member function `IndicateAtSync()`, and when `ResumeFromSync()` is called
// from charm, `AtSyncIndicator` simply passes control back to the Main chare
// via `execute_next_phase()`.
template <typename Metavariables>
class AtSyncIndicator : public CBase_AtSyncIndicator<Metavariables> {
 public:
  AtSyncIndicator(CProxy_Main<Metavariables> main_proxy);
  AtSyncIndicator(const AtSyncIndicator&) = default;
  AtSyncIndicator& operator=(const AtSyncIndicator&) = default;
  AtSyncIndicator(AtSyncIndicator&&) = default;
  AtSyncIndicator& operator=(AtSyncIndicator&&) = default;
  ~AtSyncIndicator() override {
    (void)Parallel::charmxx::RegisterChare<
        AtSyncIndicator<Metavariables>,
        CkIndex_AtSyncIndicator<Metavariables>>::registrar;
  }

  void IndicateAtSync();

  void ResumeFromSync() override;

  explicit AtSyncIndicator(CkMigrateMessage* msg)
      : CBase_AtSyncIndicator<Metavariables>(msg) {}

  void pup(PUP::er& p) override { p | main_proxy_; }

 private:
  CProxy_Main<Metavariables> main_proxy_;
};

template <typename Metavariables>
AtSyncIndicator<Metavariables>::AtSyncIndicator(
    CProxy_Main<Metavariables> main_proxy)
    : main_proxy_{main_proxy} {
  this->usesAtSync = true;
  this->setMigratable(false);
}

template <typename Metavariables>
void AtSyncIndicator<Metavariables>::IndicateAtSync() {
  this->AtSync();
}

template <typename Metavariables>
void AtSyncIndicator<Metavariables>::ResumeFromSync() {
  main_proxy_.execute_next_phase();
}
}  // namespace detail

// ================================================================

template <typename Metavariables>
Main<Metavariables>::Main(CkArgMsg* msg) {
  Informer::print_startup_info(msg);

  /// \todo detail::register_events_to_trace();

  namespace bpo = boost::program_options;
  try {
    bpo::options_description command_line_options;
    // disable clang-format because it combines the repeated call operator
    // invocations making the code more difficult to parse.
    // clang-format off
    command_line_options.add_options()
        ("help,h", "Describe program options")
        ("check-options", "Check input file options")
        ("dump-source-tree-as", bpo::value<std::string>(),
         "If specified, then a gzip archive of the source tree is dumped "
         "with the specified name. The archive can be expanded using "
         "'tar -xzf ARCHIVE.tar.gz'")
        ("dump-paths",
         "Dump the PATH, CPATH, LD_LIBRARY_PATH, LIBRARY_PATH, and "
         "CMAKE_PREFIX_PATH at compile time.")
        ("dump-environment",
         "Dump the result of printenv at compile time.")
        ("dump-build-info",
         "Dump the contents of SpECTRE's BuildInfo.txt")
        ("dump-only",
         "Exit after dumping requested information.")
        ;
    // clang-format on

    constexpr bool has_options = tmpl::size<option_list>::value > 0;
    // Add input-file option if it makes sense
    make_overloader(
        [&command_line_options](std::true_type /*meta*/, auto mv,
                                int /*gcc_bug*/)
            -> std::void_t<decltype(
                tmpl::type_from<decltype(mv)>::input_file)> {
          // Metavariables has options and default input file name
          command_line_options.add_options()(
              "input-file",
              bpo::value<std::string>()->default_value(
                  tmpl::type_from<decltype(mv)>::input_file),
              "Input file name");
        },
        [&command_line_options](std::true_type /*meta*/, auto /*mv*/,
                                auto... /*unused*/) {
          // Metavariables has options and no default input file name
          command_line_options.add_options()(
              "input-file", bpo::value<std::string>(), "Input file name");
        },
        [](std::false_type /*meta*/, auto mv, int /*gcc_bug*/)
            -> std::void_t<decltype(
                tmpl::type_from<decltype(mv)>::input_file)> {
          // Metavariables has no options and default input file name

          // always false, but must depend on mv
          static_assert(std::is_same_v<decltype(mv), void>,
                        "Metavariables supplies input file name, "
                        "but there are no options");
          ERROR("This should have failed at compile time");
        },
        [](std::false_type /*meta*/, auto... /*unused*/) {
          // Metavariables has no options and no default input file name
        })(std::bool_constant<has_options>{}, tmpl::type_<Metavariables>{}, 0);

    bpo::command_line_parser command_line_parser(msg->argc, msg->argv);
    command_line_parser.options(command_line_options);

    const bool ignore_unrecognized_command_line_options = make_overloader(
        [](auto mv, int /*gcc_bug*/)
            -> decltype(tmpl::type_from<decltype(
                            mv)>::ignore_unrecognized_command_line_options) {
          return tmpl::type_from<decltype(
              mv)>::ignore_unrecognized_command_line_options;
        },
        [](auto /*mv*/, auto... /*meta*/) { return false; })(
        tmpl::type_<Metavariables>{}, 0);
    if (ignore_unrecognized_command_line_options) {
      // Allow unknown --options
      command_line_parser.allow_unregistered();
    } else {
      // Forbid positional parameters
      command_line_parser.positional({});
    }

    bpo::variables_map parsed_command_line_options;
    bpo::store(command_line_parser.run(), parsed_command_line_options);
    bpo::notify(parsed_command_line_options);

    Options::Parser<option_list> options(Metavariables::help);

    if (parsed_command_line_options.count("help") != 0) {
      Parallel::printf("%s\n%s", command_line_options, options.help());
      sys::exit();
    }

    if (parsed_command_line_options.count("dump-source-tree-as") != 0) {
      formaline::write_to_file(
          parsed_command_line_options["dump-source-tree-as"].as<std::string>());
      Parallel::printf("Dumping archive of source tree at link time.\n");
    }
    if (parsed_command_line_options.count("dump-paths") != 0) {
      Parallel::printf("Paths at link time were:\n%s\n",
                       formaline::get_paths());
    }
    if (parsed_command_line_options.count("dump-environment") != 0) {
      Parallel::printf("Environment variables at link time were:\n%s\n",
                       formaline::get_environment_variables());
    }
    if (parsed_command_line_options.count("dump-build-info") != 0) {
      Parallel::printf("BuildInfo.txt at link time was:\n%s\n",
                       formaline::get_build_info());
    }
    if (parsed_command_line_options.count("dump-only") != 0) {
      sys::exit();
    }

    std::string input_file;
    if (has_options) {
      if (parsed_command_line_options.count("input-file") == 0) {
        ERROR("No default input file name.  Pass --input-file.");
      }
      input_file = parsed_command_line_options["input-file"].as<std::string>();
      options.parse_file(input_file);
    } else {
      options.parse("");
    }

    if (parsed_command_line_options.count("check-options") != 0) {
      // Force all the options to be created.
      options.template apply<option_list, Metavariables>([](auto... args) {
        (void)std::initializer_list<char>{((void)args, '0')...};
      });
      if (has_options) {
        Parallel::printf("\n%s parsed successfully!\n", input_file);
      } else {
        // This is still considered successful, since it means the
        // program would have started.
        Parallel::printf("\nNo options to check!\n");
      }

      // Include a check that the checkpoint dirs are available for writing as
      // part of checking the option parsing. Doing these checks together helps
      // catch more user errors before running the executable 'for real'.
      //
      // Note we don't do this check at the beginning of the Main chare
      // constructor because we don't _always_ want to error if checkpoint dirs
      // already exist. For example, running the executable with flags like
      // `--help` or `--dump-source-tree-as` should succeed even if checkpoints
      // were previously written.
      if constexpr (Algorithm_detail::has_WriteCheckpoint_v<
                        typename Metavariables::Phase>) {
        check_future_checkpoint_dirs_available();
      }

      sys::exit();
    }

    options_ =
        options.template apply<option_list, Metavariables>([](auto... args) {
          return tuples::tagged_tuple_from_typelist<option_list>(
              std::move(args)...);
        });
    Parallel::printf("\nOption parsing completed.\n");
  } catch (const bpo::error& e) {
    ERROR(e.what());
  }

  if constexpr (Algorithm_detail::has_WriteCheckpoint_v<
                    typename Metavariables::Phase>) {
    check_future_checkpoint_dirs_available();
  }

  mutable_global_cache_proxy_ = CProxy_MutableGlobalCache<Metavariables>::ckNew(
      Parallel::create_from_options<Metavariables>(
          options_, mutable_global_cache_tags{}));

  // global_cache_proxy_ depends on mutable_global_cache_proxy_.
  CkEntryOptions mutable_global_cache_dependency;
  mutable_global_cache_dependency.setGroupDepID(
      mutable_global_cache_proxy_.ckGetGroupID());

  global_cache_proxy_ = CProxy_GlobalCache<Metavariables>::ckNew(
      Parallel::create_from_options<Metavariables>(options_,
                                                   const_global_cache_tags{}),
      mutable_global_cache_proxy_, this->thisProxy,
      &mutable_global_cache_dependency);

  if constexpr (Algorithm_detail::has_LoadBalancing_v<
                    typename Metavariables::Phase>) {
    at_sync_indicator_proxy_ =
        detail::CProxy_AtSyncIndicator<Metavariables>::ckNew();
    at_sync_indicator_proxy_[0].insert(this->thisProxy, sys::my_proc());
    at_sync_indicator_proxy_.doneInserting();
  }

  tuples::tagged_tuple_from_typelist<parallel_component_tag_list>
      the_parallel_components;

  // Print info on DataBox variants
#ifdef SPECTRE_DEBUG
  Parallel::printf("\nParallel components:\n");
  tmpl::for_each<component_list>([](auto parallel_component_v) {
    using parallel_component = tmpl::type_from<decltype(parallel_component_v)>;
    using chare_type = typename parallel_component::chare_type;
    using charm_type = Parallel::charm_types_with_parameters<
        parallel_component, typename Parallel::get_array_index<
                                chare_type>::template f<parallel_component>>;
    using algorithm = typename charm_type::algorithm;
    using databox_types = typename algorithm::databox_types;
    Parallel::printf("  %s (%s) has %u DataBox variants with [",
                     pretty_type::name<parallel_component>(),
                     pretty_type::name<chare_type>(),
                     tmpl::size<databox_types>::value);
    tmpl::for_each<databox_types>([&parallel_component_v](auto databox_type_v) {
      // gcc11 and clang12 don't agree on where parallel_component_v is used...
      (void)(parallel_component_v);
      using databox_type = tmpl::type_from<decltype(databox_type_v)>;
      Parallel::printf("%u",
                       tmpl::size<typename databox_type::tags_list>::value);
      // gcc11.1.1 couldn't handle inlining the following type alias
      using last_databox_type = tmpl::back<databox_types>;
      if constexpr (not std::is_same_v<databox_type, last_databox_type>) {
        Parallel::printf(", ");
      }
    });
    Parallel::printf("] items.\n");
  });
  Parallel::printf("\n");
#endif  // SPECTRE_DEBUG

  // Construct the group proxies with a dependency on the GlobalCache proxy
  using group_component_list = tmpl::filter<
      component_list,
      tmpl::or_<Parallel::is_group_proxy<tmpl::bind<
                    Parallel::proxy_from_parallel_component, tmpl::_1>>,
                Parallel::is_node_group_proxy<tmpl::bind<
                    Parallel::proxy_from_parallel_component, tmpl::_1>>>>;
  CkEntryOptions global_cache_dependency;
  global_cache_dependency.setGroupDepID(global_cache_proxy_.ckGetGroupID());

  tmpl::for_each<group_component_list>([this, &the_parallel_components,
                                        &global_cache_dependency](
                                           auto parallel_component_v) {
    using parallel_component = tmpl::type_from<decltype(parallel_component_v)>;
    using ParallelComponentProxy =
        Parallel::proxy_from_parallel_component<parallel_component>;
    tuples::get<tmpl::type_<ParallelComponentProxy>>(the_parallel_components) =
        ParallelComponentProxy::ckNew(
            global_cache_proxy_,
            Parallel::create_from_options<Metavariables>(
                options_, typename parallel_component::initialization_tags{}),
            &global_cache_dependency);
  });

  // Create proxies for empty array chares (whose elements will be created by
  // the allocate functions of the array components during
  // execute_initialization_phase)
  using array_component_list = tmpl::filter<
      component_list,
      tmpl::and_<Parallel::is_array_proxy<tmpl::bind<
                     Parallel::proxy_from_parallel_component, tmpl::_1>>,
                 tmpl::not_<Parallel::is_bound_array<tmpl::_1>>>>;
  tmpl::for_each<array_component_list>([&the_parallel_components](
                                           auto parallel_component) {
    using ParallelComponentProxy = Parallel::proxy_from_parallel_component<
        tmpl::type_from<decltype(parallel_component)>>;
    tuples::get<tmpl::type_<ParallelComponentProxy>>(the_parallel_components) =
        ParallelComponentProxy::ckNew();
  });

  // Create proxies for empty bound array chares
  using bound_array_component_list = tmpl::filter<
      component_list,
      tmpl::and_<Parallel::is_array_proxy<tmpl::bind<
                     Parallel::proxy_from_parallel_component, tmpl::_1>>,
                 Parallel::is_bound_array<tmpl::_1>>>;
  tmpl::for_each<bound_array_component_list>([&the_parallel_components](
                                                 auto parallel_component) {
    using ParallelComponentProxy = Parallel::proxy_from_parallel_component<
        tmpl::type_from<decltype(parallel_component)>>;
    CkArrayOptions opts;
    opts.bindTo(
        tuples::get<tmpl::type_<Parallel::proxy_from_parallel_component<
            typename tmpl::type_from<decltype(parallel_component)>::bind_to>>>(
            the_parallel_components));
    tuples::get<tmpl::type_<ParallelComponentProxy>>(the_parallel_components) =
        ParallelComponentProxy::ckNew(opts);
  });

  // Send the complete list of parallel_components to the GlobalCache on
  // each Charm++ node.  After all nodes have finished, the callback is
  // executed.
  CkCallback callback(
      CkIndex_Main<Metavariables>::
          allocate_array_components_and_execute_initialization_phase(),
      this->thisProxy);
  global_cache_proxy_.set_parallel_components(the_parallel_components,
                                              callback);

  PhaseControl::initialize_phase_change_decision_data(
      make_not_null(&phase_change_decision_data_),
      *Parallel::local_branch(global_cache_proxy_));
}

template <typename Metavariables>
Main<Metavariables>::Main(CkMigrateMessage* msg)
    : CBase_Main<Metavariables>(msg) {}

template <typename Metavariables>
void Main<Metavariables>::pup(PUP::er& p) {  // NOLINT
  p | current_phase_;
  p | mutable_global_cache_proxy_;
  p | global_cache_proxy_;
  p | at_sync_indicator_proxy_;
  // Note: we do NOT serialize the options.
  // This is because options are only used in the initialization phase when
  // the executable first starts up. Thereafter, the information from the
  // options will be held in various code objects that will themselves be
  // serialized.
  p | phase_change_decision_data_;

  p | checkpoint_dir_counter_;
  if constexpr (Algorithm_detail::has_WriteCheckpoint_v<
                    typename Metavariables::Phase>) {
    if (p.isUnpacking()) {
      check_future_checkpoint_dirs_available();
    }
  }

  // For now we only support restarts on the same hardware configuration (same
  // number of nodes and same procs per node) used when writing the checkpoint.
  // We check this by adding counters to the pup stream.
  if (p.isUnpacking()) {
    int previous_nodes = 0;
    int previous_procs = 0;
    p | previous_nodes;
    p | previous_procs;
    if (previous_nodes != sys::number_of_nodes() or
        previous_procs != sys::number_of_procs()) {
      ERROR(
          "Must restart on the same hardware configuration used when writing "
          "the checkpoint.\n"
          "Checkpoint written with "
          << previous_nodes << " nodes, " << previous_procs
          << " procs.\n"
             "Restarted with "
          << sys::number_of_nodes() << " nodes, " << sys::number_of_procs()
          << " procs.");
    }
  } else {
    int current_nodes = sys::number_of_nodes();
    int current_procs = sys::number_of_procs();
    p | current_nodes;
    p | current_procs;
  }
}

template <typename Metavariables>
std::unordered_set<size_t> Main<Metavariables>::allocate_singleton_components(
    const bool avoid_proc_0) {
  ASSERT(current_phase_ == Metavariables::Phase::Initialization,
         "Must be in the Initialization phase.");
  using singleton_component_list = tmpl::filter<
      component_list,
      tmpl::bind<Parallel::is_chare_type,
                 tmpl::pin<Parallel::Algorithms::Singleton>, tmpl::_1>>;

  const size_t number_of_procs = static_cast<size_t>(sys::number_of_procs());
  const auto increment_proc = [&number_of_procs](size_t this_proc) {
    return this_proc + 1 == number_of_procs ? 0 : this_proc + 1;
  };

  std::unordered_set<size_t> procs_to_ignore{};

  size_t singleton_which_proc = 0;

  // If we should avoid proc 0, add it to the procs_to_ignore and increment the
  // starting proc number for other singletons
  if (avoid_proc_0) {
    procs_to_ignore.insert(0);
    ++singleton_which_proc;
  }

  tmpl::for_each<singleton_component_list>(
      [this, &procs_to_ignore, &singleton_which_proc,
       &increment_proc](auto singleton_component_v) {
        using singleton_component =
            tmpl::type_from<decltype(singleton_component_v)>;
        auto& local_cache = *Parallel::local_branch(global_cache_proxy_);
        auto& array_proxy =
            Parallel::get_parallel_component<singleton_component>(local_cache);
        auto options = Parallel::create_from_options<Metavariables>(
            options_, typename singleton_component::initialization_tags{});
        // If there is resource info for this singleton, use it to allocate the
        // singleton
        if constexpr (tmpl::list_contains_v<option_list,
                                            Parallel::OptionTags::SingletonInfo<
                                                singleton_component>>) {
          const auto& singleton_info =
              get<Parallel::Tags::SingletonInfo<singleton_component>>(options);
          int proc = singleton_info.proc();
          const bool exclusive = singleton_info.exclusive();
          // -1 is a sentinal value indicating that the proc should be chosen
          // automatically
          if (proc == -1) {
            while (procs_to_ignore.find(singleton_which_proc) !=
                   procs_to_ignore.end()) {
              singleton_which_proc = increment_proc(singleton_which_proc);
            }
            proc = singleton_which_proc;
          }
          // Reserve this proc so no DG elements are placed on this proc
          if (exclusive) {
            procs_to_ignore.insert(static_cast<size_t>(proc));
          }
          array_proxy[0].insert(global_cache_proxy_, std::move(options), proc);
        } else {
          // otherwise just place singletons on different procs, even if they
          // are cheap
          while (procs_to_ignore.find(singleton_which_proc) !=
                 procs_to_ignore.end()) {
            singleton_which_proc = increment_proc(singleton_which_proc);
          }
          array_proxy[0].insert(global_cache_proxy_, std::move(options),
                                singleton_which_proc);
        }
        array_proxy.doneInserting();
        singleton_which_proc = increment_proc(singleton_which_proc);
      });

  return procs_to_ignore;
}

template <typename Metavariables>
void Main<Metavariables>::
    allocate_array_components_and_execute_initialization_phase() {
  ASSERT(current_phase_ == Metavariables::Phase::Initialization,
         "Must be in the Initialization phase.");
  using array_component_list = tmpl::filter<
      component_list,
      tmpl::bind<Parallel::is_chare_type,
                 tmpl::pin<Parallel::Algorithms::Array>, tmpl::_1>>;

  bool avoid_proc_0 = false;
  // Check options to see if a component added the Parallel::Tags::AvoidProc0
  // tag to its initialization tags. If one did, get the value.
  if constexpr (tmpl::list_contains_v<option_list,
                                      Parallel::OptionTags::AvoidProc0>) {
    const auto avoid_proc_0_input =
        Parallel::create_from_options<Metavariables>(
            options_, tmpl::list<Parallel::Tags::AvoidProc0>{});
    avoid_proc_0 = get<Parallel::Tags::AvoidProc0>(avoid_proc_0_input);
  }

  // Allocate all singleton components, specifying if proc 0 should be avoided
  std::unordered_set<size_t> procs_to_ignore =
      allocate_singleton_components(avoid_proc_0);

  tmpl::for_each<array_component_list>([this, &procs_to_ignore](
                                           auto parallel_component_v) {
    using parallel_component = tmpl::type_from<decltype(parallel_component_v)>;
    // This is a Spectre array component built on a Charm++ array chare. The
    // component is in charge of allocating and distributing its elements over
    // the computing system.
    auto local_options = Parallel::create_from_options<Metavariables>(
        options_, typename parallel_component::initialization_tags{});
    parallel_component::allocate_array(global_cache_proxy_, local_options,
                                       procs_to_ignore);
  });

  // Free any resources from the initial option parsing.
  options_ = decltype(options_){};

  tmpl::for_each<component_list>([this](auto parallel_component_v) {
    using parallel_component = tmpl::type_from<decltype(parallel_component_v)>;
    Parallel::get_parallel_component<parallel_component>(
        *Parallel::local_branch(global_cache_proxy_))
        .start_phase(current_phase_);
  });
  CkStartQD(CkCallback(CkIndex_Main<Metavariables>::execute_next_phase(),
                       this->thisProxy));
}

template <typename Metavariables>
void Main<Metavariables>::execute_next_phase() {
  current_phase_ = Metavariables::determine_next_phase(
      make_not_null(&phase_change_decision_data_), current_phase_,
      global_cache_proxy_);
  if (Metavariables::Phase::Exit == current_phase_) {
    Informer::print_exit_info();
    sys::exit();
  }
  tmpl::for_each<component_list>([this](auto parallel_component) {
    tmpl::type_from<decltype(parallel_component)>::execute_next_phase(
        current_phase_, global_cache_proxy_);
  });

  // Here we handle phases with direct Charm++ calls. By handling these phases
  // after calling each component's execute_next_phase entry method, we ensure
  // that each component knows what phase it is in. This is useful for pup
  // functions that need special handling that depends on the phase.
  //
  // Note that in future versions of Charm++ it may become possible for pup
  // functions to have knowledge of the migration type. At that point, it
  // should no longer be necessary to wait until after
  // component::execute_next_phase to make the direct charm calls. Instead, the
  // load balance or checkpoint work could be initiated *before* the call to
  // component::execute_next_phase and *without* the need for a quiescence
  // detection. This may be a slight optimization.
  if constexpr (Algorithm_detail::has_LoadBalancing_v<
                    typename Metavariables::Phase>) {
    if (current_phase_ == Metavariables::Phase::LoadBalancing) {
      CkStartQD(CkCallback(CkIndex_Main<Metavariables>::start_load_balance(),
                           this->thisProxy));
      return;
    }
  }
  if constexpr (Algorithm_detail::has_WriteCheckpoint_v<
                    typename Metavariables::Phase>) {
    if (current_phase_ == Metavariables::Phase::WriteCheckpoint) {
      CkStartQD(
          CkCallback(CkIndex_Main<Metavariables>::start_write_checkpoint(),
                     this->thisProxy));
      return;
    }
  }

  // The general case simply returns to execute_next_phase
  CkStartQD(CkCallback(CkIndex_Main<Metavariables>::execute_next_phase(),
                       this->thisProxy));
}

template <typename Metavariables>
void Main<Metavariables>::start_load_balance() {
  at_sync_indicator_proxy_.IndicateAtSync();
  // No need for a callback to return to execute_next_phase: this is done by
  // ResumeFromSync instead.
}

template <typename Metavariables>
void Main<Metavariables>::start_write_checkpoint() {
  const std::string dir = checkpoint_dir();
  checkpoint_dir_counter_++;
  CkStartCheckpoint(
      dir.c_str(), CkCallback(CkIndex_Main<Metavariables>::execute_next_phase(),
                              this->thisProxy));
}

template <typename Metavariables>
template <typename InvokeCombine, typename... Tags>
void Main<Metavariables>::phase_change_reduction(
    ReductionData<ReductionDatum<tuples::TaggedTuple<Tags...>, InvokeCombine,
                                 funcl::Identity, std::index_sequence<>>>
        reduction_data) {
  using tagged_tuple_type = std::decay_t<
      std::tuple_element_t<0, std::decay_t<decltype(reduction_data.data())>>>;
  (void)Parallel::charmxx::RegisterPhaseChangeReduction<
      Metavariables, InvokeCombine, Tags...>::registrar;
  static_assert(tt::is_a_v<tuples::TaggedTuple, tagged_tuple_type>,
                "The main chare expects a tagged tuple in the phase change "
                "reduction target.");
  reduction_data.finalize();
  PhaseControl::TaggedTupleMainCombine::apply(
      make_not_null(&phase_change_decision_data_),
      get<0>(reduction_data.data()));
}

/// @{
/// Send data from a parallel component to the Main chare for making
/// phase-change decisions.
///
/// This function is distinct from `Parallel::contribute_to_reduction` because
/// this function contributes reduction data to the Main chare via the entry
/// method `phase_change_reduction`. This must be done because the entry method
/// must alter member data specific to the Main chare, and the Main chare cannot
/// execute actions like most SpECTRE parallel components.
/// For all cases other than sending phase-change decision data to the Main
/// chare, you should use `Parallel::contribute_to_reduction`.
template <typename SenderComponent, typename ArrayIndex, typename Metavariables,
          class... Ts>
void contribute_to_phase_change_reduction(
    tuples::TaggedTuple<Ts...> data_for_reduction,
    Parallel::GlobalCache<Metavariables>& cache,
    const ArrayIndex& array_index) {
  using phase_change_tags_and_combines_list =
      PhaseControl::get_phase_change_tags<Metavariables>;
  using reduction_data_type =
      PhaseControl::reduction_data<tmpl::list<Ts...>,
                                   phase_change_tags_and_combines_list>;
  (void)Parallel::charmxx::RegisterReducerFunction<
      reduction_data_type::combine>::registrar;
  CkCallback callback(CProxy_Main<Metavariables>::index_t::
                          template redn_wrapper_phase_change_reduction<
                              PhaseControl::TaggedTupleCombine, Ts...>(nullptr),
                      cache.get_main_proxy().value());
  reduction_data_type reduction_data{data_for_reduction};
  Parallel::local(
      Parallel::get_parallel_component<SenderComponent>(cache)[array_index])
      ->contribute(static_cast<int>(reduction_data.size()),
                   reduction_data.packed().get(),
                   Parallel::charmxx::charm_reducer_functions.at(
                       std::hash<Parallel::charmxx::ReducerFunctions>{}(
                           &reduction_data_type::combine)),
                   callback);
}
template <typename SenderComponent, typename Metavariables, class... Ts>
void contribute_to_phase_change_reduction(
    tuples::TaggedTuple<Ts...> data_for_reduction,
    Parallel::GlobalCache<Metavariables>& cache) {
  using phase_change_tags_and_combines_list =
      PhaseControl::get_phase_change_tags<Metavariables>;
  using reduction_data_type =
      PhaseControl::reduction_data<tmpl::list<Ts...>,
                                   phase_change_tags_and_combines_list>;
  (void)Parallel::charmxx::RegisterReducerFunction<
      reduction_data_type::combine>::registrar;
  CkCallback callback(CProxy_Main<Metavariables>::index_t::
                          template redn_wrapper_phase_change_reduction<
                              PhaseControl::TaggedTupleCombine, Ts...>(nullptr),
                      cache.get_main_proxy().value());
  reduction_data_type reduction_data{data_for_reduction};
  // Note that Singletons could be supported by directly calling the main
  // entry function, but due to this and other peculiarities with
  // Singletons, it is best to discourage their use.
  static_assert(
      not std::is_same_v<typename SenderComponent::chare_type,
                         Parallel::Algorithms::Singleton>,
      "Phase change reduction is not supported for singleton chares. "
      "Consider constructing your chare as a length-1 array chare if you "
      "need to contribute to phase change data");
  Parallel::local_branch(
      Parallel::get_parallel_component<SenderComponent>(cache))
      ->contribute(static_cast<int>(reduction_data.size()),
                   reduction_data.packed().get(),
                   Parallel::charmxx::charm_reducer_functions.at(
                       std::hash<Parallel::charmxx::ReducerFunctions>{}(
                           &reduction_data_type::combine)),
                   callback);
}
/// @}

template <typename Metavariables>
std::tuple<std::string, std::string, size_t>
Main<Metavariables>::checkpoint_dir_basename_pad() const {
  const std::string basename = "SpectreCheckpoint";
  constexpr size_t pad = 6;

  const std::string counter = std::to_string(checkpoint_dir_counter_);
  const std::string padded_counter =
      std::string(pad - counter.size(), '0').append(counter);
  const std::string checkpoint_dir = basename + padded_counter;
  return std::make_tuple(checkpoint_dir, basename, pad);
}

template <typename Metavariables>
std::string Main<Metavariables>::checkpoint_dir() const {
  const auto [checkpoint_dir, basename, pad] = checkpoint_dir_basename_pad();
  (void)basename;
  (void)pad;
  if (file_system::check_if_dir_exists(checkpoint_dir)) {
    ERROR("Can't write checkpoint: dir " + checkpoint_dir + " already exists!");
  }
  return checkpoint_dir;
}

template <typename Metavariables>
void Main<Metavariables>::check_future_checkpoint_dirs_available() const {
  // Can't lambda-capture from structured binding in clang-11
  std::string checkpoint_dir;
  std::string basename;
  size_t pad;
  std::tie(checkpoint_dir, basename, pad) = checkpoint_dir_basename_pad();

  // Find existing files with names that match the checkpoint dir name pattern
  const auto all_files = file_system::ls();
  const std::regex re(basename + "[0-9]{" + std::to_string(pad) + "}");
  std::vector<std::string> checkpoint_files;
  std::copy_if(all_files.begin(), all_files.end(),
               std::back_inserter(checkpoint_files),
               [&re](const std::string& s) { return std::regex_match(s, re); });

  // Using string comparison of filenames, check that all the files we found
  // are from older checkpoints, but not from future checkpoints
  const bool found_older_checkpoints_only = std::all_of(
      checkpoint_files.begin(), checkpoint_files.end(),
      [&checkpoint_dir](const std::string& s) { return s < checkpoint_dir; });
  if (not found_older_checkpoints_only) {
    ERROR(
        "Can't start run: found checkpoints that may be overwritten!\n"
        "Dirs from "
        << checkpoint_dir << " onward must not exist.\n");
  }
}

}  // namespace Parallel

#define CK_TEMPLATES_ONLY
#include "Parallel/Main.def.h"
#undef CK_TEMPLATES_ONLY
