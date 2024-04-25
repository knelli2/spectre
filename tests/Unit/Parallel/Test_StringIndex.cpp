// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "Parallel/StringIndex.hpp"

SPECTRE_TEST_CASE("Unit.Parallel.StringIndex", "[Unit][Parallel]") {
  const std::string name_1{"Name1"};
  const std::string name_2{"Name2"};

  const StringIndex index_1{name_1};
  const StringIndex index_2{name_2};

  CHECK(index_1 == index_1);
  CHECK_FALSE(index_1 == index_2);
  CHECK_FALSE(index_1 != index_1);
  CHECK(index_1 != index_2);

  StringIndex fill_1{};
  StringIndex fill_2{};

  std::memset(&fill_1, 0, sizeof(fill_1));
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 7)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif  // defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 7)
  std::memset(&fill_2, 255, sizeof(fill_2));
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 7)
#pragma GCC diagnostic pop
#endif  // defined(__GNUC__) && !defined(__clang__) && (__GNUC__ > 7)

  new (&fill_1) StringIndex(name_1);
  new (&fill_2) StringIndex(name_2);

  using Hash = std::hash<StringIndex>;

  CHECK((fill_1 == fill_2) == (index_1 == index_2));
  CHECK((fill_1 != fill_2) == (index_1 != index_2));
  CHECK((Hash{}(fill_1) == Hash{}(fill_2)) == (index_1 == index_2));
}
