// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Framework/TestCreation.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Parallel/Tags/ResourceInfo.hpp"

namespace Parallel {
namespace {
struct ParallelComponent {};
}  // namespace

SPECTRE_TEST_CASE("Unit.Parallel.Tags.ResourceInfo", "[Unit][Parallel]") {
  TestHelpers::db::test_simple_tag<Tags::SingletonInfo<ParallelComponent>>(
      "ParallelComponent");
  TestHelpers::db::test_simple_tag<Tags::AvoidProc0>("AvoidProc0");

  // Has to be 0, otherwise an error will occur because each test is run on one
  // core and there is a check in SingletonInfo that the input proc isn't beyond
  // the "last" proc
  const auto singleton_info = TestHelpers::test_option_tag<
      OptionTags::SingletonInfo<ParallelComponent>>(
      "Proc: 0\n"
      "Exclusive: false\n");
  const auto singleton_info_auto = TestHelpers::test_option_tag<
      OptionTags::SingletonInfo<ParallelComponent>>(
      "Proc: Auto\n"
      "Exclusive: true\n");

  CHECK(singleton_info.proc() == 0);
  CHECK_FALSE(singleton_info.exclusive());
  // -1 is a sentinal value for auto
  CHECK(singleton_info_auto.proc() == -1);
  CHECK(singleton_info_auto.exclusive());

  const bool avoid_proc_0_f =
      TestHelpers::test_option_tag<OptionTags::AvoidProc0>("false");
  const bool avoid_proc_0_t =
      TestHelpers::test_option_tag<OptionTags::AvoidProc0>("true");

  CHECK(avoid_proc_0_t);
  CHECK_FALSE(avoid_proc_0_f);
}

// [[OutputRegex, non-negative integer]]
SPECTRE_TEST_CASE("Unit.Parallel.Tags.ResourceInfo.BadNegativeProc",
                  "[Unit][Parallel]") {
  ERROR_TEST();
  [[maybe_unused]] const auto singleton_info = TestHelpers::test_option_tag<
      OptionTags::SingletonInfo<ParallelComponent>>(
      "Proc: -2\n"
      "Exclusive: false\n");
}

// [[OutputRegex, is beyond the last proc]]
SPECTRE_TEST_CASE("Unit.Parallel.Tags.ResourceInfo.BadProcBeyondLast",
                  "[Unit][Parallel]") {
  ERROR_TEST();
  [[maybe_unused]] const auto singleton_info = TestHelpers::test_option_tag<
      OptionTags::SingletonInfo<ParallelComponent>>(
      "Proc: 2\n"
      "Exclusive: false\n");
}

// [[OutputRegex, singleton has requested to be placed on proc 0]]
SPECTRE_TEST_CASE("Unit.Parallel.Tags.ResourceInfo.BadInconsistentProc0",
                  "[Unit][Parallel]") {
  ERROR_TEST();
  const Parallel::SingletonInfo<ParallelComponent> singleton_info{{0}, false};
  const bool avoid_proc_0 = true;

  [[maybe_unused]] const auto created_singleton_info =
      Parallel::Tags::SingletonInfo<ParallelComponent>::create_from_options(
          singleton_info, avoid_proc_0);
}
}  // namespace Parallel
