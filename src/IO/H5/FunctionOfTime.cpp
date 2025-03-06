// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "IO/H5/FunctionOfTime.hpp"

#include <algorithm>
#include <hdf5.h>
#include <iosfwd>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "IO/H5/CheckH5.hpp"
#include "IO/H5/Header.hpp"
#include "IO/H5/Helpers.hpp"
#include "IO/H5/OpenGroup.hpp"
#include "IO/H5/Type.hpp"
#include "IO/H5/Version.hpp"
#include "IO/H5/Wrappers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/Serialize.hpp"
#include "Utilities/StdHelpers.hpp"

namespace h5::detail {
FunctionOfTime::FunctionOfTime(
    const bool exists, detail::OpenGroup&& group, const hid_t location,
    const std::string& name,
    const std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>&
        initial_function_of_time,
    const uint32_t version)
    : group_(std::move(group)),
      name_(extension() == name.substr(name.size() > extension().size()
                                           ? name.size() - extension().size()
                                           : 0)
                ? name
                : name + extension()),
      path_(group_.group_path_with_trailing_slash() + name),
      version_(version) {
  if (exists) {
    (void)location;
    ERROR("Cannot support reading function of time from H5 yet.");
    {
      // We treat this as an internal version for now. We'll need to deal with
      // proper versioning later.

      // Check if the version exists before calling the open_version.
      const htri_t version_exists = H5Aexists(dataset_id_, "version.ver");
      if (version_exists != 0) {
        const Version open_version(true, detail::OpenGroup{}, dataset_id_,
                                   "version");
        version_ = open_version.get_version();
      }
    }
    {
      const htri_t header_exists = H5Aexists(dataset_id_, "header.hdr");
      if (header_exists != 0) {
        const Header header(true, detail::OpenGroup{}, dataset_id_, "header");
        header_ = header.get_header();
      }
    }
  } else {  // file does not exist
    if (initial_function_of_time == nullptr) {
      ERROR(
          "Must supply an initial function of time when creating a new .fot "
          "subfile");
    }

    const auto serialized_initial_function_of_time =
        serialize(initial_function_of_time);

    // Write initial FoT
    write_data(group_.id(), serialized_initial_function_of_time,
               {serialized_initial_function_of_time.size()},
               "initial_function_of_time"s, true);

    // Write time bounds
    // QUESTION: Does std::numeric_limits<double>::infinity() work here?
    const std::array<double, 2> time_bounds =
        initial_function_of_time->time_bounds();
    write_data(group_.id(), std::vector{time_bounds[0], time_bounds[1]},
               std::vector{2_st}, "time_bounds"s, true);

    {
      Version open_version(false, detail::OpenGroup{}, dataset_id_, "version",
                           version_);
    }
    {
      Header header(false, detail::OpenGroup{}, dataset_id_, "header");
      header_ = header.get_header();
    }
  }
}

void FunctionOfTime::append(const double time, const DataVector& update) {
  if (not updates_.has_value()) {
    // TODO: What should this be?
    hid_t location = -1;
    // QUESTION: New group??
    std::vector<std::string> legend(update.size() + 1);
    legend[0] = "Time";
    for (size_t i = 0; i < update.size(); i++) {
      legend[i + 1] = "Component_" + std::to_string(i);
    }
    updates_ =
        h5::Dat{false, group_, location, "Updates", std::move(legend), 0};
  } else if (const size_t expected_legend_size =
                 updates_.value().get_legend().size();
             expected_legend_size != update.size() + 1) {
    ERROR("Update for H5 function of time is incorrect size. Expecting "
          << expected_legend_size << " but got " << update.size() + 1);
  }
}
}  // namespace h5::detail
