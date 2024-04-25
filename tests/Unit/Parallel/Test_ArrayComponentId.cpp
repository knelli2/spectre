// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <ckarrayindex.h>
#include <functional>
#include <type_traits>

#include "Domain/Structure/ElementId.hpp"
#include "Framework/TestHelpers.hpp"
#include "Parallel/ArrayComponentId.hpp"
#include "Parallel/ArrayIndex.hpp"
#include "Parallel/StringIndex.hpp"

namespace {
struct Component0 {};
struct Component1 {};

SPECTRE_TEST_CASE("Unit.Parallel.ArrayComponentId", "[Unit][Parallel]") {
  using ArrayComponentId = Parallel::ArrayComponentId;
  using Hash = std::hash<ArrayComponentId>;
  ArrayComponentId c0id0{std::add_pointer_t<Component0>{nullptr},
                         Parallel::ArrayIndex<ElementId<1>>(ElementId<1>{0})};
  ArrayComponentId c0id1{std::add_pointer_t<Component0>{nullptr},
                         Parallel::ArrayIndex<ElementId<1>>(ElementId<1>{1})};
  const std::string string_index_name{"KylesBallerIndex"};
  ArrayComponentId c0str{
      std::add_pointer_t<Component0>{nullptr},
      Parallel::ArrayIndex<StringIndex>{StringIndex{string_index_name}}};

  ArrayComponentId c1id0 =
      Parallel::make_array_component_id<Component1>(ElementId<1>{0});
  ArrayComponentId c1id1 =
      Parallel::make_array_component_id<Component1>(ElementId<1>{1});
  ArrayComponentId c1str =
      Parallel::make_array_component_id<Component1>(string_index_name);

  CHECK(c0id0 == c0id0);
  CHECK_FALSE(c0id0 != c0id0);
  CHECK(c0id0 != c0id1);
  CHECK(c0id0 != c0str);
  CHECK(c0id0 != c1id0);
  CHECK(c0id0 != c1id1);
  CHECK(c0id0 != c1str);

  CHECK(c0id0.component_id() == c0id0.component_id());
  CHECK_FALSE(c0id0.component_id() != c0id0.component_id());
  CHECK(c0id0.component_id() == c0id1.component_id());
  CHECK(c0id0.component_id() == c0str.component_id());
  CHECK(c0id0.component_id() != c1id0.component_id());
  CHECK(c0id0.component_id() != c1id1.component_id());
  CHECK(c0id0.component_id() != c1str.component_id());
  CHECK(c0str.component_id() != c1str.component_id());

  CHECK(c0id0.array_index() == c0id0.array_index());
  CHECK_FALSE(c0id0.array_index() == c0id1.array_index());
  CHECK_FALSE(c0id0.array_index() == c0str.array_index());
  CHECK(c0id0.array_index() == c1id0.array_index());
  CHECK_FALSE(c0id0.array_index() == c1id1.array_index());
  CHECK_FALSE(c0id0.array_index() == c1str.array_index());
  CHECK(c0str.array_index() == c1str.array_index());

  CHECK(Hash{}(c0id0) == Hash{}(c0id0));
  CHECK(Hash{}(c0id0) != Hash{}(c0id1));
  CHECK(Hash{}(c0id0) != Hash{}(c1id0));
  CHECK(Hash{}(c0id0) != Hash{}(c1id1));
  CHECK(Hash{}(c0id0) != Hash{}(c0str));
  CHECK(Hash{}(c0id0) != Hash{}(c1str));
  CHECK(Hash{}(c1id0) == Hash{}(c1id0));
  CHECK(Hash{}(c1id0) != Hash{}(c1id1));
  CHECK(Hash{}(c0str) == Hash{}(c0str));
  CHECK(Hash{}(c0str) != Hash{}(c1str));
  CHECK(Hash{}(c1str) == Hash{}(c1str));

  // Test PUP
  test_serialization(c0id0);
  test_serialization(c0id1);
  test_serialization(c0str);
  test_serialization(c1id0);
  test_serialization(c1id1);
  test_serialization(c1str);
}
}  // namespace
