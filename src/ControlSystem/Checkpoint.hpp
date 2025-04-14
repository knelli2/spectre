// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <exception>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Tags/SystemTags.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/Object.hpp"
#include "IO/H5/OpenGroup.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace control_system {
class Checkpoint : public h5::Object {
 public:
  /// \cond HIDDEN_SYMBOLS
  static std::string extension() { return ".cscp"; }

  Checkpoint(bool exists, h5::detail::OpenGroup&& group, hid_t location,
             const std::string& name, uint32_t version = 1);

  Checkpoint(const Checkpoint& /*rhs*/) = delete;
  Checkpoint& operator=(const Checkpoint& /*rhs*/) = delete;
  Checkpoint(Checkpoint&& /*rhs*/) = delete;             // NOLINT
  Checkpoint& operator=(Checkpoint&& /*rhs*/) = delete;  // NOLINT

  ~Checkpoint() override = default;
  /// \endcond HIDDEN_SYMBOLS

  void write_checkpoint(const std::vector<char>& averager,
                        const std::vector<char>& tuner,
                        const std::vector<char>& controller,
                        const std::vector<char>& control_error_class) const;

  std::unordered_map<std::string, std::vector<char>> read_checkpoint() const;

  /*!
   * \returns the header of the Checkpoint file
   */
  const std::string& get_header() const { return header_; }

  /*!
   * \returns the user-specified version number of the Checkpoint file
   *
   * \note h5::Version returns a uint32_t, so we return one here too for the
   * version
   */
  uint32_t get_version() const { return version_; }

  const std::string& subfile_path() const override { return path_; }

 private:
  /// \cond HIDDEN_SYMBOLS
  h5::detail::OpenGroup group_;
  std::string name_;
  std::string path_;
  uint32_t version_;
  h5::detail::OpenGroup serialized_group_;
  std::string header_;
  /// \endcond HIDDEN_SYMBOLS
};

template <size_t AveragerDerivOrder, bool AllowDecrease,
          size_t ControllerDerivOrder, typename ControlError>
void write_checkpoint(const Checkpoint& checkpoint_subfile,
                      const Averager<AveragerDerivOrder>& averager,
                      const TimescaleTuner<AllowDecrease>& tuner,
                      const Controller<ControllerDerivOrder>& controller,
                      const ControlError& control_error_class) {
  std::unordered_map<std::string, std::vector<char>> serialized_objects{};
  const std::vector<char> serialized_averager = serialize(averager);
  const std::vector<char> serialized_tuner = serialize(tuner);
  const std::vector<char> serialized_controller = serialize(controller);
  const std::vector<char> serialized_control_error_class =
      serialize(control_error_class);

  checkpoint_subfile.write_checkpoint(serialized_averager, serialized_tuner,
                                      serialized_controller,
                                      serialized_control_error_class);
}

template <size_t AveragerDerivOrder, bool AllowDecrease,
          size_t ControllerDerivOrder, typename ControlError>
void read_checkpoint(
    const Checkpoint& checkpoint_subfile,
    const gsl::not_null<Averager<AveragerDerivOrder>*> averager,
    const gsl::not_null<TimescaleTuner<AllowDecrease>*> tuner,
    const gsl::not_null<Controller<ControllerDerivOrder>*> controller,
    const gsl::not_null<ControlError*> control_error_class) {
  const auto serialized_classes = checkpoint_subfile.read_checkpoint();

  const auto deserialize_class = [&serialized_classes]<typename T>(
                                     const gsl::not_null<T*>& read_into,
                                     const std::string& name) {
    try {
      (*read_into) = deserialize<T>(serialized_classes.at(name).data());
    } catch (const std::exception& /*exception*/) {
      ERROR("Unable to read control system checkpoint for" << name);
    }
  };

  deserialize_class(averager, "Averager");
  deserialize_class(tuner, "Tuner");
  deserialize_class(controller, "Controller");
  deserialize_class(control_error_class, "ControlErrorClass");
}
}  // namespace control_system
