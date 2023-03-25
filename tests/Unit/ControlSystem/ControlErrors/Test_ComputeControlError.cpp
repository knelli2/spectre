// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <optional>
#include <string>

#include "ControlSystem/ControlErrors/ComputeControlError.hpp"
#include "ControlSystem/Protocols/ControlError.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/Tags.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Helpers/ControlSystem/SystemHelpers.hpp"
#include "Helpers/ControlSystem/TestStructs.hpp"
#include "Options/Options.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace control_system {
namespace {
struct SystemNoArgs : tt::ConformsTo<protocols::ControlSystem> {
  static constexpr size_t deriv_order = 2;
  static std::string name() { return "SystemNoArgs"; }
  static std::optional<std::string> component_name(
      const size_t i, const size_t /*num_components*/) {
    return get_output(i);
  }
  using measurement =
      TestHelpers::Measurement<TestHelpers::TestStructs_detail::LabelA>;
  using simple_tags = tmpl::list<>;
  using control_error = TestHelpers::ControlError<1>;
  struct process_measurement {
    using argument_tags = tmpl::list<>;
  };
};

struct TagInt : db::SimpleTag {
  using type = int;
};

struct TagDouble : db::SimpleTag {
  using type = double;
};

struct TagString : db::SimpleTag {
  using type = std::string;
};

struct ControlErrorArgs : tt::ConformsTo<protocols::ControlError> {
  static constexpr size_t expected_number_of_excisions = 0;
  using object_centers = domain::object_list<>;
  void pup(PUP::er& p) { p | set_int; }

  using options = tmpl::list<>;
  static constexpr Options::String help{""};

  ControlErrorArgs() = default;
  ControlErrorArgs(const int in) : set_int(in) {}

  using return_tags = tmpl::list<TagInt, TagDouble>;
  using argument_tags = tmpl::list<TagString>;

  template <typename Metavariables, typename... QueueTags>
  DataVector operator()(
      const gsl::not_null<int*> mutable_int,
      const gsl::not_null<double*> mutable_double,
      const std::string& arg_string,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const double /*time*/, const std::string& /*function_of_time_name*/,
      const tuples::TaggedTuple<QueueTags...>& /*measurements*/) {
    *mutable_int = static_cast<int>(arg_string.size());
    set_int = *mutable_int;
    *mutable_double = 1.234;
    return DataVector{*mutable_double};
  }

  int set_int;
};

struct SystemArgs : tt::ConformsTo<protocols::ControlSystem> {
  static constexpr size_t deriv_order = 2;
  static std::string name() { return "SystemArgs"; }
  static std::optional<std::string> component_name(
      const size_t i, const size_t /*num_components*/) {
    return get_output(i);
  }
  using measurement =
      TestHelpers::Measurement<TestHelpers::TestStructs_detail::LabelA>;
  using simple_tags = tmpl::list<>;
  using control_error = ControlErrorArgs;
  struct process_measurement {
    using argument_tags = tmpl::list<>;
  };
};

struct Metavariables {
  using component_list = tmpl::list<>;
};

void test() {
  // Avoid unused variable warnings for the constexpr vars
  (void)SystemNoArgs::deriv_order;
  (void)ControlErrorArgs::expected_number_of_excisions;
  (void)ControlErrorArgs::help;
  (void)SystemArgs::deriv_order;

  Parallel::GlobalCache<Metavariables> cache{};
  auto box = db::create<db::AddSimpleTags<TagInt, TagDouble, TagString,
                                          Tags::ControlError<SystemNoArgs>,
                                          Tags::ControlError<SystemArgs>>>(
      -5, 0.0, "HelloThere"s, SystemNoArgs::control_error{},
      SystemArgs::control_error{0});
  tuples::TaggedTuple<> empty_tuple{};

  const double time = 5.7;  // Not used anywhere but needed to pass in

  const DataVector control_error_no_args = db::mutate_apply<ComputeControlError<
      SystemNoArgs, typename SystemNoArgs::control_error::return_tags,
      typename SystemNoArgs::control_error::argument_tags>>(
      make_not_null(&box), cache, time, pretty_type::name<SystemNoArgs>(),
      empty_tuple);

  CHECK(control_error_no_args == DataVector{});
  CHECK(db::get<TagInt>(box) == -5);
  CHECK(db::get<TagDouble>(box) == 0.0);
  CHECK(db::get<TagString>(box) == "HelloThere");
  CHECK(db::get<Tags::ControlError<SystemArgs>>(box).set_int == 0);

  const DataVector control_error_args = db::mutate_apply<ComputeControlError<
      SystemArgs, typename SystemArgs::control_error::return_tags,
      typename SystemArgs::control_error::argument_tags>>(
      make_not_null(&box), cache, time, pretty_type::name<SystemArgs>(),
      empty_tuple);

  CHECK(control_error_args == DataVector{1.234});
  CHECK(db::get<TagInt>(box) == 10);
  CHECK(db::get<TagDouble>(box) == 1.234);
  CHECK(db::get<TagString>(box) == "HelloThere");
  CHECK(db::get<Tags::ControlError<SystemArgs>>(box).set_int == 10);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.ControlSystem.ControlErrors.ComputeControlError",
                  "[ControlSystem][Unit]") {
  test();
}
}  // namespace control_system
