// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Access.hpp"
#include "DataStructures/DataBox/AsAccess.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Framework/TestHelpers.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/AccessWrapper.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Callback.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Target.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/NoSuchType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct DoubleTag : db::SimpleTag {
  using type = double;
};

// Note that we are only testing the compile time interface of this class so we
// don't inherit from the Callback base class here
struct TestCallback : tt::ConformsTo<intrp2::protocols::Callback> {
  using tags_to_observe_on_target = tmpl::list<DoubleTag>;
  using non_observation_tags_on_target = tmpl::list<>;
  using volume_compute_tags = tmpl::list<>;

  using options = tmpl::list<>;
};

// Just for the DataBox compute below. Needs this type alias in the Points
// struct.
struct FakePoints {
  using tags_on_target = tmpl::list<>;
};

struct TestTarget : public tt::ConformsTo<intrp2::protocols::Target> {
  using temporal_id_tag = DoubleTag;
  using frame = NoSuchType;
  using points = FakePoints;

  using possible_runtime_callbacks = tmpl::list<TestCallback>;
  using compile_time_callbacks = tmpl::list<>;
};

using Derived = intrp2::TargetAccessWrapper<TestTarget, 3>;

struct TestMetavars {
  struct factory_creation {
    using factory_classes =
        tmpl::map<tmpl::pair<intrp2::AccessWrapper, tmpl::list<Derived>>>;
  };
};

void test() {
  register_factory_classes_with_charm<TestMetavars>();

  std::unique_ptr<intrp2::AccessWrapper> access_wrapper =
      std::make_unique<Derived>();

  const double expected_double = 9876.54321;

  db::mutate<DoubleTag>(
      [&expected_double](const gsl::not_null<double*> double_value) {
        (*double_value) = expected_double;
      },
      make_not_null(&access_wrapper->access()));

  auto serialized_access_wrapper = serialize_and_deserialize(access_wrapper);

  CHECK(db::get<DoubleTag>(serialized_access_wrapper->access()) ==
        expected_double);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ParallelAlgorithms.Interpolation.Runtime.AccessWrapper",
                  "[Unit]") {
  test();
}
