// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ControlSystem/Checkpoint.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include "IO/H5/AccessType.hpp"
#include "IO/H5/Header.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/OpenGroup.hpp"
#include "IO/H5/Version.hpp"

namespace control_system {
Checkpoint::Checkpoint(const bool exists, h5::detail::OpenGroup&& group,
                       const hid_t /*location*/, const std::string& name,
                       const uint32_t version)
    : group_(std::move(group)),
      name_(extension() == name.substr(name.size() > extension().size()
                                           ? name.size() - extension().size()
                                           : 0)
                ? name
                : name + extension()),
      path_(group.group_path_with_trailing_slash() + name),
      version_(version),
      serialized_group_(group_.id(), name_, h5::AccessType::ReadWrite) {
  if (exists) {
    {
      // We treat this as an internal version for now. We'll need to deal with
      // proper versioning later.

      // Check if the version exists before calling the open_version.
      const htri_t version_exists =
          H5Aexists(serialized_group_.id(), "version.ver");
      if (version_exists != 0) {
        const h5::Version open_version(true, h5 ::detail::OpenGroup{},
                                       serialized_group_.id(), "version");
        version_ = open_version.get_version();
      }
    }
    {
      const htri_t header_exists =
          H5Aexists(serialized_group_.id(), "header.hdr");
      if (header_exists != 0) {
        const h5::Header header(true, h5::detail::OpenGroup{},
                                serialized_group_.id(), "header");
        header_ = header.get_header();
      }
    }
  } else {
    {
      h5::Version open_version(false, h5::detail::OpenGroup{},
                               serialized_group_.id(), "version", version_);
    }
    {
      h5::Header header(false, h5::detail::OpenGroup{}, serialized_group_.id(),
                        "header");
      header_ = header.get_header();
    }
  }
}

void Checkpoint::write_checkpoint(
    const std::vector<char>& averager, const std::vector<char>& tuner,
    const std::vector<char>& controller,
    const std::vector<char>& control_error_class) const {
  h5::write_data(serialized_group_.id(), averager, {averager.size()},
                 "Averager");
  h5::write_data(serialized_group_.id(), tuner, {tuner.size()}, "Tuner");
  h5::write_data(serialized_group_.id(), controller, {controller.size()},
                 "Controller");
  // Only write the control error if there's actually data
  if (control_error_class.size() > 0) {
    h5::write_data(serialized_group_.id(), control_error_class,
                   {control_error_class.size()}, "ControlErrorClass");
  }
}

std::unordered_map<std::string, std::vector<char>> Checkpoint::read_checkpoint()
    const {
  const auto read_serialized_class = [this](const std::string& name) {
    if (not h5::contains_dataset_or_group(serialized_group_.id(), "", name)) {
      ERROR("Cannot read control system checkpoint. Missing " << name);
    }

    return h5::read_data<1, std::vector<char>>(serialized_group_.id(), name);
  };

  std::unordered_map<std::string, std::vector<char>> result{};

  for (const std::string& name : {"Averager", "Tuner", "Controller"}) {
    result[name] = read_serialized_class(name);
  }

  // This may not exist if the control error class doesn't store any data
  if (h5::contains_dataset_or_group(serialized_group_.id(), "",
                                    "ControlErrorClass")) {
    result["ControlErrorClass"] = read_serialized_class("ControlErrorClass");
  }

  return result;
}
}  // namespace control_system
