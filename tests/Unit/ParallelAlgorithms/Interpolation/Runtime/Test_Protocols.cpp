// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/NoSuchType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct TagForTarget : db::SimpleTag {
  using type = std::optional<double>;
};

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

struct TestPoints : tt::ConformsTo<intrp2::protocols::Points> {
  using tags_on_target = tmpl::list<TagForTarget>;
  using points_volume_compute_tags = tmpl::list<DataVectorComputeTag<100>>;

  using options = tmpl::list<>;
};

struct TestTarget;

// Note that we are only testing the compile time compliance of this class so we
// don't inherit from the Callback base class here
template <size_t Index>
struct TestCallback : tt::ConformsTo<intrp2::protocols::Callback> {
  using tags_to_observe_on_target = tmpl::list<DoubleTag<1>, DoubleTag<2>>;
  using non_observation_tags_on_target = tmpl::list<DataVectorTag<1>>;
  using volume_compute_tags = tmpl::list<DataVectorComputeTag<2>>;

  using options = tmpl::list<>;
};

struct TestTarget : public tt::ConformsTo<intrp2::protocols::Target> {
  using temporal_id_tag = DoubleTag<0>;
  using frame = NoSuchType;
  using points = TestPoints;

  using possible_runtime_callbacks =
      tmpl::list<TestCallback<0>, TestCallback<1>>;
  using compile_time_callbacks = tmpl::list<TestCallback<2>>;
};
}  // namespace

// Check points protocol conformance
static_assert(tt::assert_conforms_to_v<TestPoints, intrp2::protocols::Points>);

// Check callback protocol conformance
static_assert(
    tt::assert_conforms_to_v<TestCallback<0>, intrp2::protocols::Callback>);

// Check target protocol conformance
static_assert(tt::assert_conforms_to_v<TestTarget, intrp2::protocols::Target>);
