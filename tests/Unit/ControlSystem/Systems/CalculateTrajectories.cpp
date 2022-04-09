// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "IO/H5/AccessType.hpp"
#include "IO/H5/Dat.hpp"
#include "IO/H5/File.hpp"
#include "IO/H5/Object.hpp"
#include "Utilities/FileSystem.hpp"
#include "Utilities/Gsl.hpp"

namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame

namespace {
using ExpansionMap = domain::CoordinateMaps::TimeDependent::CubicScale<3>;
using RotationMap = domain::CoordinateMaps::TimeDependent::Rotation<3>;

using CoordMap = domain::CoordinateMap<Frame::Grid, Frame::Inertial,
                                       ExpansionMap, RotationMap>;

SPECTRE_TEST_CASE("Unit.ApparentHorizon.Trajectories",
                  "[Unit][ApparentHorizons]") {
  domain::FunctionsOfTime::register_derived_with_charm();
  const std::string filename{"/home/knelli/GhBinaryBlackHoleReductionData.h5"};

  const h5::H5File<h5::AccessType::ReadOnly> my_file(filename);
  const auto& rotation_x_dat =
      my_file.get<h5::Dat>("/ControlSystems/Rotation/X");
  const auto& rotation_y_dat =
      my_file.get<h5::Dat>("/ControlSystems/Rotation/Y");
  const auto& rotation_z_dat =
      my_file.get<h5::Dat>("/ControlSystems/Rotation/Z");
  const auto& expansion_dat = my_file.get<h5::Dat>("/ControlSystems/Expansion");
  const Matrix rotation_x_data = rotation_x_dat.get_data();
  const Matrix rotation_y_data = rotation_y_dat.get_data();
  const Matrix rotation_z_data = rotation_z_dat.get_data();
  const Matrix expansion_data = expansion_dat.get_data();

  const double init_time = 0.0;
  std::array<DataVector, 3> init_expansion{
      {DataVector{1, 1.0}, DataVector{1, -4.6148457646200002e-05},
       DataVector{1, 0.0}}};
  std::array<DataVector, 3> init_rotation{
      DataVector{3, 0.0}, DataVector{{0.0, 0.0, 1.5264577062000000e-02}},
      DataVector{3, 0.0}};
  std::array<DataVector, 1> init_quat{{DataVector{{1.0, 0.0, 0.0, 0.0}}}};

  domain::FunctionsOfTime::PiecewisePolynomial<2> exp_fot{
      init_time, init_expansion, expansion_data(0, 0)};
  domain::FunctionsOfTime::QuaternionFunctionOfTime<2> rot_fot{
      init_time, init_quat, init_rotation, rotation_z_data(0, 0)};

  const size_t num_updates =
      std::min({expansion_data.rows(), rotation_x_data.rows(),
                rotation_y_data.rows(), rotation_z_data.rows()});

  for (size_t i = 0; i < num_updates - 1; i++) {
    const double exp_update_time = expansion_data(i, 0);
    const double rot_update_time = rotation_x_data(i, 0);
    const double exp_expr_time = expansion_data(i + 1, 0);
    const double rot_expr_time = rotation_x_data(i + 1, 0);

    const DataVector exp_updated_deriv{{expansion_data(i, 3)}};
    const DataVector rot_updated_deriv{
        {rotation_x_data(i, 3), rotation_y_data(i, 3), rotation_x_data(i, 3)}};

    exp_fot.update(exp_update_time, exp_updated_deriv, exp_expr_time);
    rot_fot.update(rot_update_time, rot_updated_deriv, rot_expr_time);
  }

  const std::string expansion_name{"Expansion"};
  const std::string rotation_name{"Rotation"};

  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      functions_of_time{};
  functions_of_time[expansion_name] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          exp_fot);
  // values taken from input file
  functions_of_time[expansion_name + "OuterBoundary"] =
      std::make_unique<domain::FunctionsOfTime::FixedSpeedCubic>(
          1.0, init_time, -1.0e-6, -1.0e-6);
  functions_of_time[rotation_name] =
      std::make_unique<domain::FunctionsOfTime::QuaternionFunctionOfTime<2>>(
          rot_fot);

  ExpansionMap expansion_map{300.0, expansion_name,
                             expansion_name + "OuterBoundary"s};
  RotationMap rotation_map{rotation_name};

  CoordMap coord_map{expansion_map, rotation_map};

  const double x_pos = 7.683;
  tnsr::I<double, 3, Frame::Grid> init_pos_A{};
  tnsr::I<double, 3, Frame::Grid> init_pos_B{};
  get<0>(init_pos_A) = -x_pos;
  get<0>(init_pos_B) = x_pos;
  for (size_t i = 1; i < 3; i++) {
    init_pos_A.get(i) = 0.0;
    init_pos_B.get(i) = 0.0;
  }

  const std::vector<std::string> legend{
      {"Time", "InertialCenter_X", "InertialCenter_Y", "InertialCenter_Z"}};

  const std::string outfile_name{"/home/knelli/trajectories.h5"};
  if (file_system::check_if_file_exists(outfile_name)) {
    file_system::rm(outfile_name, true);
  }
  h5::H5File<h5::AccessType::ReadWrite> out_file(outfile_name, true);
  auto& traj_AhA =
      out_file.template insert<h5::Dat>("/HorizonTrajectories/AhA", legend);
  auto& traj_AhB =
      out_file.template insert<h5::Dat>("/HorizonTrajectories/AhB", legend);

  for (size_t i = 0; i < num_updates; i++) {
    const double time = rotation_x_data(i, 0);

    const auto pos_AhA = coord_map(init_pos_A, time, functions_of_time);
    const auto pos_AhB = coord_map(init_pos_B, time, functions_of_time);

    std::vector<double> vec_pos_AhA{4, 0.0};
    std::vector<double> vec_pos_AhB{4, 0.0};

    vec_pos_AhA[0] = time;
    vec_pos_AhB[0] = time;

    for (size_t j = 0; j < 3; j++) {
      vec_pos_AhA[j + 1] = pos_AhA.get(j);
      vec_pos_AhB[j + 1] = pos_AhB.get(j);
    }

    traj_AhA.append(vec_pos_AhA);
    traj_AhB.append(vec_pos_AhB);
  }
}
}  // namespace
