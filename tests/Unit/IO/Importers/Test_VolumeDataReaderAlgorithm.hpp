// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#define CATCH_CONFIG_RUNNER

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/VolumeData.hpp"
#include "IO/Importers/Actions/ReadVolumeData.hpp"
#include "IO/Importers/Actions/ReceiveVolumeData.hpp"
#include "IO/Importers/Actions/RegisterWithElementDataReader.hpp"
#include "IO/Importers/ElementDataReader.hpp"
#include "IO/Importers/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/Algorithms/AlgorithmArray.hpp"
#include "Parallel/Algorithms/AlgorithmSingleton.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Local.hpp"
#include "Parallel/Main.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "ParallelAlgorithms/Actions/SetData.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

struct ScalarFieldTag : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "ScalarField"; }
};

template <size_t Dim>
struct VectorFieldTag : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim>;
  static std::string name() { return "VectorField"; }
};

enum class Grid { Fine, Coarse };

std::ostream& operator<<(std::ostream& os, Grid grid) {
  switch (grid) {
    case Grid::Fine:
      return os << "Fine";
    case Grid::Coarse:
      return os << "Coarse";
    default:
      ERROR("Missing case for grid");
  }
}

template <Grid TheGrid>
constexpr size_t number_of_elements = 2;
template <>
constexpr size_t number_of_elements<Grid::Coarse> = 1;

/// [option_group]
template <Grid TheGrid>
struct VolumeDataOptions {
  using group = importers::OptionTags::Group;
  static std::string name() { return MakeString{} << TheGrid << "VolumeData"; }
  static constexpr Options::String help = "Numeric volume data";
};
/// [option_group]

template <size_t Dim, Grid TheGrid>
struct TestVolumeData {
  static constexpr size_t num_elements = number_of_elements<TheGrid>;
  std::array<std::array<size_t, Dim>, num_elements> extents;
  std::array<DataVector, num_elements> scalar_field_data;
  std::array<std::array<DataVector, Dim>, num_elements> vector_field_data;
};

template <size_t Dim>
const TestVolumeData<Dim, Grid::Fine> fine_volume_data{};
template <>
const TestVolumeData<1, Grid::Fine> fine_volume_data<1>{
    {{// Grid extents on first element
      {{2}},
      // Grid extents on second element
      {{3}}}},
    {{// Field on first element
      {{1., 2.}},
      // Field on second element
      {{3., 4., 5.}}}},
    {{// Vector components on first element
      {{{{1., 2.}}}},
      // Vector components on second element
      {{{{7., 8., 9.}}}}}}};
template <>
const TestVolumeData<2, Grid::Fine> fine_volume_data<2>{
    {{{{2, 2}}, {{3, 2}}}},
    {{{{1., 2., 3., 4.}}, {{3., 4., 5., 6., 7., 8.}}}},
    {{{{{{1., 2., 3., 4.}}, {{3., 4., 5., 6.}}}},
      {{{{7., 8., 9., 10., 11., 12.}}, {{10., 11., 12., 13., 14., 15.}}}}}}};
template <>
const TestVolumeData<3, Grid::Fine> fine_volume_data<3>{
    {{{{2, 2, 2}}, {{3, 2, 2}}}},
    {{{{1., 2., 3., 4., 5., 6., 7., 8.}},
      {{3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.}}}},
    {{{{{{1., 2., 3., 4., 5., 6., 7., 8.}},
        {{3., 4., 5., 6., 7., 8., 9., 10.}},
        {{5., 6., 7., 8., 9., 10., 11., 12.}}}},
      {{{{7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18.}},
        {{10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21.}},
        {{13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.}}}}}}};

template <size_t Dim>
const TestVolumeData<Dim, Grid::Coarse> coarse_volume_data{};
template <>
const TestVolumeData<1, Grid::Coarse> coarse_volume_data<1>{
    {{{{2}}}}, {{{{17., 18.}}}}, {{{{{{19., 20.}}}}}}};
template <>
const TestVolumeData<2, Grid::Coarse> coarse_volume_data<2>{
    {{{{2, 2}}}},
    {{{{17., 18., 19., 20.}}}},
    {{{{{{19., 20., 21., 22.}}, {{21., 22., 23., 24.}}}}}}};
template <>
const TestVolumeData<3, Grid::Coarse> coarse_volume_data<3>{
    {{{{2, 2, 2}}}},
    {{{{17., 18., 19., 20., 21., 22., 23., 24.}}}},
    {{{{{{19., 20., 21., 22., 23., 24., 25., 26.}},
        {{21., 22., 23., 24., 25., 26., 27., 28.}},
        {{23., 24., 25., 26., 27., 28., 29., 30.}}}}}}};

template <bool Check>
void clean_test_data(const std::string& data_file_name) {
  if (file_system::check_if_file_exists(data_file_name)) {
    file_system::rm(data_file_name, true);
  } else if (Check) {
    ERROR("Expected test data file '" << data_file_name << "' does not exist");
  }
}

template <size_t Dim, Grid TheGrid>
void write_test_data(const std::string& data_file_name,
                     const std::string& subgroup,
                     const double observation_value,
                     const TestVolumeData<Dim, TheGrid>& test_data) {
  // Open file for test data
  h5::H5File<h5::AccessType::ReadWrite> data_file{data_file_name, true};
  auto& test_data_file = data_file.insert<h5::VolumeData>("/" + subgroup);

  // Construct test data for all elements
  std::vector<ElementVolumeData> element_data{};
  for (size_t i = 0; i < number_of_elements<TheGrid>; i++) {
    const std::string element_name = MakeString{} << ElementId<Dim>{i};
    std::vector<TensorComponent> tensor_components{};
    std::vector<size_t> element_extents{};
    std::vector<Spectral::Quadrature> element_quadratures{};
    std::vector<Spectral::Basis> element_bases{};
    tensor_components.push_back(
        {element_name + "/ScalarField", test_data.scalar_field_data[i]});
    for (size_t d = 0; d < Dim; d++) {
      static const std::array<std::string, 3> dim_suffix{{"x", "y", "z"}};
      tensor_components.push_back(
          {element_name + "/VectorField_" + dim_suffix[d],
           test_data.vector_field_data[i][d]});
      element_extents.push_back(test_data.extents[i][d]);
      element_quadratures.push_back(Spectral::Quadrature::Gauss);
      element_bases.push_back(Spectral::Basis::Chebyshev);
    }
    element_data.push_back(
        {std::move(element_extents), std::move(tensor_components),
         std::move(element_bases), std::move(element_quadratures)});
  }
  test_data_file.write_volume_data(0, observation_value,
                                   std::move(element_data));
}

template <size_t Dim>
struct WriteTestData {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    // The data may be in a shared file, so first clean both, then write both
    clean_test_data<false>(
        get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Fine>>>(box));
    clean_test_data<false>(
        get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Coarse>>>(box));
    write_test_data<Dim>(
        get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Fine>>>(box),
        get<importers::Tags::Subgroup<VolumeDataOptions<Grid::Fine>>>(box),
        std::get<double>(get<importers::Tags::ObservationValue<
                             VolumeDataOptions<Grid::Fine>>>(box)),
        fine_volume_data<Dim>);
    write_test_data<Dim>(
        get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Coarse>>>(box),
        get<importers::Tags::Subgroup<VolumeDataOptions<Grid::Coarse>>>(box),
        std::get<double>(get<importers::Tags::ObservationValue<
                             VolumeDataOptions<Grid::Coarse>>>(box)),
        coarse_volume_data<Dim>);
    return {std::move(box), true};
  }
};

struct CleanTestData {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    clean_test_data<true>(
        get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Fine>>>(box));
    if (get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Coarse>>>(box) !=
        get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Fine>>>(box)) {
      clean_test_data<true>(
          get<importers::Tags::FileGlob<VolumeDataOptions<Grid::Coarse>>>(box));
    }
    return {std::move(box), true};
  }
};

template <size_t Dim, typename Metavariables>
struct TestDataWriter {
  using chare_type = Parallel::Algorithms::Singleton;
  using metavariables = Metavariables;
  using phase_dependent_action_list =
      tmpl::list<Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::Initialization,
                                        tmpl::list<WriteTestData<Dim>>>,

                 Parallel::PhaseActions<typename Metavariables::Phase,
                                        Metavariables::Phase::TestResult,
                                        tmpl::list<CleanTestData>>>;
  using initialization_tags = tmpl::list<>;

  static void initialize(
      Parallel::CProxy_GlobalCache<Metavariables>& /*global_cache*/) {}

  static void execute_next_phase(
      const typename Metavariables::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_component = Parallel::get_parallel_component<TestDataWriter>(
        *Parallel::local_branch(global_cache));
    local_component.start_phase(next_phase);
  }
};

template <size_t Dim>
struct InitializeElement {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ElementId<Dim>& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) {
    return std::make_tuple(
        ::Initialization::merge_into_databox<
            InitializeElement,
            db::AddSimpleTags<ScalarFieldTag, VectorFieldTag<Dim>>>(
            std::move(box), Scalar<DataVector>{}, tnsr::I<DataVector, Dim>{}),
        true);
  }
};

template <size_t Dim, Grid TheGrid>
void test_result(const ElementId<Dim>& element_index,
                 const TestVolumeData<Dim, TheGrid>& test_data,
                 const Scalar<DataVector>& scalar_field,
                 const tnsr::I<DataVector, Dim>& vector_field) {
  const size_t raw_element_index = ElementId<Dim>{element_index}.block_id();
  Scalar<DataVector> expected_scalar_field{
      test_data.scalar_field_data[raw_element_index]};
  SPECTRE_PARALLEL_REQUIRE(scalar_field == expected_scalar_field);
  tnsr::I<DataVector, Dim> expected_vector_field{};
  for (size_t d = 0; d < Dim; d++) {
    expected_vector_field[d] =
        test_data.vector_field_data[raw_element_index][d];
  }
  SPECTRE_PARALLEL_REQUIRE(vector_field == expected_vector_field);
}

template <size_t Dim, Grid TheGrid>
struct TestResult {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    if (TheGrid == Grid::Fine) {
      test_result(array_index, fine_volume_data<Dim>, get<ScalarFieldTag>(box),
                  get<VectorFieldTag<Dim>>(box));
    } else {
      test_result(array_index, coarse_volume_data<Dim>,
                  get<ScalarFieldTag>(box), get<VectorFieldTag<Dim>>(box));
    }
    return {std::move(box), true};
  }
};

template <size_t Dim, Grid TheGrid, typename Metavariables>
struct ElementArray {
  using chare_type = Parallel::Algorithms::Array;
  using array_index = ElementId<Dim>;
  using metavariables = Metavariables;
  using initialization_tags = tmpl::list<>;
  using array_allocation_tags = tmpl::list<>;

  using import_fields = tmpl::list<ScalarFieldTag, VectorFieldTag<Dim>>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::Initialization,
                             tmpl::list<InitializeElement<Dim>>>,
      /// [import_actions]
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Register,
          tmpl::list<importers::Actions::RegisterWithElementDataReader,
                     Parallel::Actions::TerminatePhase>>,
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::ImportData,
          tmpl::list<importers::Actions::ReadVolumeData<
                         VolumeDataOptions<TheGrid>, import_fields>,
                     importers::Actions::ReceiveVolumeData<
                         VolumeDataOptions<TheGrid>, import_fields>,
                     Parallel::Actions::TerminatePhase>>,
      /// [import_actions]
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::TestResult,
                             tmpl::list<TestResult<Dim, TheGrid>>>>;

  static void allocate_array(
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache,
      const tuples::tagged_tuple_from_typelist<initialization_tags>&
          initialization_items) {
    auto& array_proxy = Parallel::get_parallel_component<ElementArray>(
        *Parallel::local_branch(global_cache));

    for (size_t i = 0, which_proc = 0,
                number_of_procs = static_cast<size_t>(sys::number_of_procs());
         i < number_of_elements<TheGrid>; i++) {
      ElementId<Dim> element_index{i};
      array_proxy[element_index].insert(global_cache, initialization_items,
                                        which_proc);
      which_proc = which_proc + 1 == number_of_procs ? 0 : which_proc + 1;
    }
    array_proxy.doneInserting();
  }

  static void execute_next_phase(
      const typename Metavariables::Phase next_phase,
      Parallel::CProxy_GlobalCache<Metavariables>& global_cache) {
    auto& local_cache = *Parallel::local_branch(global_cache);
    Parallel::get_parallel_component<ElementArray>(local_cache)
        .start_phase(next_phase);
  }
};

/// [metavars]
template <size_t Dim>
struct Metavariables {
  using component_list =
      tmpl::list<ElementArray<Dim, Grid::Fine, Metavariables>,
                 ElementArray<Dim, Grid::Coarse, Metavariables>,
                 TestDataWriter<Dim, Metavariables>,
                 importers::ElementDataReader<Metavariables>>;

  static constexpr const char* const help{"Test the volume data reader"};
  static constexpr bool ignore_unrecognized_command_line_options = false;

  enum class Phase { Initialization, Register, ImportData, TestResult, Exit };

  template <typename... Tags>
  static Phase determine_next_phase(
      const gsl::not_null<
          tuples::TaggedTuple<Tags...>*> /*phase_change_decision_data*/,
      const Phase& current_phase,
      const Parallel::CProxy_GlobalCache<Metavariables>& /*cache_proxy*/) {
    switch (current_phase) {
      case Phase::Initialization:
        return Phase::Register;
      case Phase::Register:
        return Phase::ImportData;
      case Phase::ImportData:
        return Phase::TestResult;
      default:
        return Phase::Exit;
    }
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) {}
};
/// [metavars]

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling, &setup_memory_allocation_failure_reporting};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
