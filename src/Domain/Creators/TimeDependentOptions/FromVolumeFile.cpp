// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/TimeDependentOptions/FromVolumeFile.hpp"

#include <array>
#include <boost/math/quaternion.hpp>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/QuaternionHelpers.hpp"
#include "Domain/Structure/ObjectLabel.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/VolumeData.hpp"
#include "Options/Context.hpp"
#include "Options/ParseError.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace domain::creators::time_dependent_options {

namespace {
// Returns a clone of the FoT requested
std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> get_function_of_time(
    const std::string& function_of_time_name, const std::string& h5_filename,
    const std::string& subfile_name, const double time,
    const Options::Context& context) {
  const h5::H5File<h5::AccessType::ReadOnly> h5_file{h5_filename};
  const auto& vol_file = h5_file.get<h5::VolumeData>(subfile_name);

  const std::vector<size_t> obs_ids = vol_file.list_observation_ids();
  if (obs_ids.empty()) {
    PARSE_ERROR(context, function_of_time_name
                             << ": There are no observation IDs in the subfile "
                             << subfile_name << " of H5 file " << h5_filename);
  }
  // Take last observation ID so we have all possible times available
  std::optional<std::vector<char>> serialized_functions_of_time =
      vol_file.get_functions_of_time(obs_ids[obs_ids.size() - 1]);

  if (not serialized_functions_of_time.has_value()) {
    PARSE_ERROR(context,
                function_of_time_name
                    << ": There are no functions of time in the subfile "
                    << subfile_name << " of the H5 file " << h5_filename
                    << ". Choose a different subfile or H5 file.");
  }

  const auto functions_of_time = deserialize<std::unordered_map<
      std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>>(
      serialized_functions_of_time->data());

  if (not functions_of_time.contains(function_of_time_name)) {
    PARSE_ERROR(context, "No function of time named "
                             << function_of_time_name << " in the subfile "
                             << subfile_name << " of the H5 file "
                             << h5_filename);
  }

  const auto& function_of_time = functions_of_time.at(function_of_time_name);
  const std::array<double, 2> time_bounds = function_of_time->time_bounds();

  if (time < time_bounds[0] or time > time_bounds[1]) {
    using ::operator<<;
    PARSE_ERROR(context, function_of_time_name
                             << ": The requested time " << time
                             << " is out of the range of the function of time "
                             << time_bounds << " from the subfile "
                             << subfile_name << " of the H5 file "
                             << h5_filename);
  }

  return function_of_time->get_clone();
}
}  // namespace

template <typename FoTName>
FromVolumeFile<FoTName>::FromVolumeFile(const std::string& h5_filename,
                                        const std::string& subfile_name,
                                        const double time,
                                        const Options::Context& context) {
  const std::string function_of_time_name = pretty_type::name<FoTName>();

  const auto function_of_time = get_function_of_time(
      function_of_time_name, h5_filename, subfile_name, time, context);

  values = function_of_time->func_and_2_derivs(time);
}

FromVolumeFile<names::Expansion>::FromVolumeFile(
    const std::string& h5_filename, const std::string& subfile_name,
    const double time, const Options::Context& context) {
  const std::string expansion_function_of_time_name{"Expansion"};

  const auto exp_function_of_time =
      get_function_of_time(expansion_function_of_time_name, h5_filename,
                           subfile_name, time, context);
  const auto exp_outer_boundary_function_of_time =
      get_function_of_time(expansion_function_of_time_name + "OuterBoundary",
                           h5_filename, subfile_name, time, context);

  expansion_values = exp_function_of_time->func_and_2_derivs(time);
  expansion_values_outer_boundary =
      exp_outer_boundary_function_of_time->func_and_2_derivs(time);
}

FromVolumeFile<names::Rotation>::FromVolumeFile(
    const std::string& h5_filename, const std::string& subfile_name,
    const double time, const bool is_quat_fot,
    const Options::Context& context) {
  const std::string function_of_time_name{"Rotation"};

  const auto function_of_time = get_function_of_time(
      function_of_time_name, h5_filename, subfile_name, time, context);

  if (is_quat_fot) {
    bool filled_values = false;
    const auto set_values = [this, &function_of_time, &time,
                             &filled_values]<size_t N>() {
      const auto* const quat_function_of_time =
          dynamic_cast<domain::FunctionsOfTime::QuaternionFunctionOfTime<N>*>(
              function_of_time.get());

      if (quat_function_of_time != nullptr) {
        quaternions = quat_function_of_time->quat_func_and_2_derivs(time);
        angle_values = quat_function_of_time->angle_func_and_2_derivs(time);
        filled_values = true;
      }
    };

    set_values.operator()<2>();
    set_values.operator()<3>();

    if (not filled_values) {
      PARSE_ERROR(context,
                  function_of_time_name
                      << ": The max deriv of the quaternion function of "
                         "time must be either 2 or 3.");
    }
  } else {
    quaternions = function_of_time->func_and_2_derivs(time);
    std::array<boost::math::quaternion<double>, 3> quats;
    for (size_t i = 0; i < quaternions.size(); i++) {
      gsl::at(quats, i) = datavector_to_quaternion(gsl::at(quaternions, i));
    }

    // We can't compute the angle so set it to zero
    angle_values[0] = DataVector{4, 0.0};
    angle_values[1] = quaternion_to_datavector(2.0 * conj(quats[0]) * quats[1]);
    angle_values[2] = quaternion_to_datavector(
        conj(quats[0]) * (2.0 * quats[2] - quats[1] * datavector_to_quaternion(
                                                          angle_values[0])));
    // We'd need more quaternion derivatives to compute the higher angle
    // derivatives. Also the above DataVectors have size 4, so we need to make
    // them have size 3.
    angle_values[2] = DataVector{4, 0.0};
    for (size_t i = 0; i < 3; i++) {
      gsl::at(angle_values, i) =
          DataVector{gsl::at(angle_values, i)[1], gsl::at(angle_values, i)[2],
                     gsl::at(angle_values, i)[3]};
    }

    // angle_values = std::nullopt;
  }
}

template <ObjectLabel Object>
FromVolumeFileShapeSize<Object>::FromVolumeFileShapeSize(
    const std::string& h5_filename, const std::string& subfile_name,
    const double time, const Options::Context& context) {
  const std::string object = domain::name(Object);
  const std::string shape_function_of_time_name{"Shape" + object};
  const std::string size_function_of_time_name{"Size" + object};

  const auto shape_function_of_time = get_function_of_time(
      shape_function_of_time_name, h5_filename, subfile_name, time, context);
  const auto size_function_of_time = get_function_of_time(
      size_function_of_time_name, h5_filename, subfile_name, time, context);

  shape_values = shape_function_of_time->func_and_2_derivs(time);
  size_values = size_function_of_time->func_and_2_derivs(time);
}

template struct FromVolumeFile<names::Translation>;
template struct FromVolumeFileShapeSize<ObjectLabel::A>;
template struct FromVolumeFileShapeSize<ObjectLabel::B>;
template struct FromVolumeFileShapeSize<ObjectLabel::None>;
}  // namespace domain::creators::time_dependent_options
