// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <optional>
#include <type_traits>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Metafunctions.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/NoSuchType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <size_t Index>
struct DoubleTag : db::SimpleTag {
  using type = double;
};
template <size_t Index>
struct DataVectorTag : db::SimpleTag {
  using type = DataVector;
};
template <size_t Index>
struct DataVectorComputeTag : DataVectorTag<Index>, db::ComputeTag {
  using base = DataVectorTag<Index>;
  using return_type = DataVector;
  using argument_tags = DoubleTag<Index>;

  static void function(const gsl::not_null<DataVector*> dv,
                       const double argument) {
    *dv = DataVector{argument};
  }
};

struct TestTarget;

// Note that we are only testing the compile time compliance of this class so we
// don't inherit from the Callback base class here
template <size_t Index>
struct TestCallback : tt::ConformsTo<intrp2::protocols::Callback> {
  using tags_to_observe_on_target =
      tmpl::list<DoubleTag<Index>, DataVectorComputeTag<Index>,
                 DataVectorComputeTag<Index + 10>>;
  using non_observation_tags_on_target = tmpl::list<DoubleTag<Index + 10>>;
  using volume_compute_tags = tmpl::list<>;

  using options = tmpl::list<>;
};

struct TestTarget : public tt::ConformsTo<intrp2::protocols::Target> {
  using temporal_id_tag = DoubleTag<0>;
  using frame = NoSuchType;
  using points = NoSuchType;

  using possible_runtime_callbacks =
      tmpl::list<TestCallback<0>, TestCallback<1>>;
  using compile_time_callbacks = tmpl::list<TestCallback<2>>;
};

template <typename List1, typename List2>
struct check_lists
    : tmpl::and_<typename std::is_same<tmpl::list_difference<List1, List2>,
                                       tmpl::list<>>::type,
                 typename std::is_same<tmpl::list_difference<List1, List2>,
                                       tmpl::list<>>::type>::type {};

using expected_simple_tags =
    tmpl::list<DoubleTag<0>, DataVectorTag<0>, DataVectorTag<10>>;
using simple_tags = intrp2::metafunctions::simple_tags_from_mixed_tags<
    typename TestCallback<0>::tags_to_observe_on_target>;
static_assert(check_lists<expected_simple_tags, simple_tags>::value);

using expected_callbacks =
    tmpl::list<TestCallback<0>, TestCallback<1>, TestCallback<2>>;
using callbacks = intrp2::metafunctions::all_callbacks<TestTarget>;
static_assert(check_lists<expected_callbacks, callbacks>::value);

using expected_tags_to_observe =
    tmpl::list<DoubleTag<0>, DataVectorTag<0>, DataVectorTag<10>, DoubleTag<1>,
               DataVectorTag<1>, DataVectorTag<11>, DoubleTag<2>,
               DataVectorTag<2>, DataVectorTag<12>>;
using tags_to_observe = intrp2::metafunctions::all_tags_to_observe<TestTarget>;
static_assert(check_lists<expected_tags_to_observe, tags_to_observe>::value);
}  // namespace
