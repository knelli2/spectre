// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <limits>
#include <memory>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Checkpoint.hpp"
#include "ControlSystem/ControlErrors/Size.hpp"
#include "ControlSystem/ControlErrors/Size/Initial.hpp"
#include "ControlSystem/ControlErrors/Size/RegisterDerivedWithCharm.hpp"
#include "ControlSystem/ControlErrors/Size/State.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/CheckH5.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/Header.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/OpenGroup.hpp"
#include "IO/H5/Version.hpp"
#include "IO/H5/Wrappers.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace {

void test() {
  const std::string filename{"Haircut.h5"};
  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }

  const Averager<2> averager{0.25, true};
  const TimescaleTuner<false> tuner{{0.1}, 10.0, 1.e-4, 1e-3, 1.01};
  const Controller<2> controller{0.3};
  const control_system::ControlErrors::Size<2, ::domain::ObjectLabel::None>
      size_control_error{
          4, 0.25,
          TimescaleTuner<true>{{0.1}, 20.0, 1.e-4, 2.5e-4, 1.01, 1.e-3, 0.99},
          std::make_unique<control_system::size::States::Initial>(),
          std::nullopt};

  const std::vector<char> expected_serialized_averager = serialize(averager);
  const std::vector<char> expected_serialized_tuner = serialize(tuner);
  const std::vector<char> expected_serialized_controller =
      serialize(controller);
  const std::vector<char> expected_serialized_size_control_error =
      serialize(size_control_error);

  {
    h5::H5File<h5::AccessType::ReadWrite> h5_file{filename};

    {
      auto& checkpoint_subfile =
          h5_file.insert<control_system::Checkpoint>("ControlSystem1");

      CHECK(checkpoint_subfile.extension() == ".cscp");
      CHECK(checkpoint_subfile.get_version() == 1);
      CHECK(checkpoint_subfile.subfile_path() == "/ControlSystem1");

      checkpoint_subfile.write_checkpoint(
          expected_serialized_averager, expected_serialized_tuner,
          expected_serialized_controller,
          expected_serialized_size_control_error);

      h5_file.close_current_object();
    }
    {
      auto& checkpoint_subfile =
          h5_file.insert<control_system::Checkpoint>("ControlSystem2");

      CHECK(checkpoint_subfile.extension() == ".cscp");
      CHECK(checkpoint_subfile.get_version() == 1);
      CHECK(checkpoint_subfile.subfile_path() == "/ControlSystem2");

      control_system::write_checkpoint(checkpoint_subfile, averager, tuner,
                                       controller, size_control_error);

      h5_file.close_current_object();
    }
  }

  {
    const h5::H5File<h5::AccessType::ReadOnly> h5_file{filename};

    REQUIRE(h5_file.exists<control_system::Checkpoint>("ControlSystem1"));
    REQUIRE(h5_file.exists<control_system::Checkpoint>("ControlSystem2"));

    {
      const auto& checkpoint_subfile =
          h5_file.get<control_system::Checkpoint>("ControlSystem2");

      const auto serialized_checkpoint = checkpoint_subfile.read_checkpoint();
      CHECK(serialized_checkpoint.at("Averager") ==
            expected_serialized_averager);
      CHECK(serialized_checkpoint.at("Tuner") == expected_serialized_tuner);
      CHECK(serialized_checkpoint.at("Controller") ==
            expected_serialized_controller);
      CHECK(serialized_checkpoint.at("ControlErrorClass") ==
            expected_serialized_size_control_error);

      h5_file.close_current_object();
    }

    {
      const auto& checkpoint_subfile =
          h5_file.get<control_system::Checkpoint>("ControlSystem1");

      Averager<2> deserialized_averager{};
      TimescaleTuner<false> deserialized_tuner{};
      Controller<2> deserialized_controller{};
      control_system::ControlErrors::Size<2, ::domain::ObjectLabel::None>
          deserialized_size_control_error{};

      control_system::read_checkpoint(
          checkpoint_subfile, make_not_null(&deserialized_averager),
          make_not_null(&deserialized_tuner),
          make_not_null(&deserialized_controller),
          make_not_null(&deserialized_size_control_error));

      CHECK(deserialized_averager == averager);
      CHECK(deserialized_tuner == tuner);
      CHECK(deserialized_controller == controller);
      // No equality operator for size control error so we call a couple
      // functions
      CHECK(
          deserialized_size_control_error.discontinuous_change_has_occurred() ==
          size_control_error.discontinuous_change_has_occurred());
      CHECK(deserialized_size_control_error.get_suggested_timescale() ==
            size_control_error.get_suggested_timescale());

      h5_file.close_current_object();
    }
  }

  if (file_system::check_if_file_exists(filename)) {
    file_system::rm(filename, true);
  }
}

SPECTRE_TEST_CASE("Unit.ControlSystem.Checkpoint", "[Unit][ControlSystem]") {
  control_system::size::register_derived_with_charm();
  test();
}
}  // namespace
