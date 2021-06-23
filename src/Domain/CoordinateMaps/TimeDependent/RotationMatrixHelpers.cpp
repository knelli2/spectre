// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/CoordinateMaps/TimeDependent/RotationMatrixHelpers.hpp"

#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace {
void AddBilinearTerm(Matrix& rot_matrix, DataVector q1, DataVector q2,
                     double coef = 1.0) {
  rot_matrix(0, 0) +=
      coef * (q1[0] * q2[0] + q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]);
  rot_matrix(1, 1) +=
      coef * (q1[0] * q2[0] + q1[2] * q2[2] - q1[1] * q2[1] - q1[3] * q2[3]);
  rot_matrix(2, 2) +=
      coef * (q1[0] * q2[0] + q1[3] * q2[3] - q1[1] * q2[1] - q1[2] * q2[2]);

  rot_matrix(0, 1) += coef * (2.0 * (q1[1] * q2[2] - q1[0] * q2[3]));
  rot_matrix(0, 2) += coef * (2.0 * (q1[0] * q2[2] + q1[1] * q2[3]));
  rot_matrix(1, 0) += coef * (2.0 * (q1[1] * q2[2] + q1[0] * q2[3]));
  rot_matrix(1, 2) += coef * (2.0 * (q1[2] * q2[3] - q1[0] * q2[1]));
  rot_matrix(2, 0) += coef * (2.0 * (q1[1] * q2[3] - q1[0] * q2[2]));
  rot_matrix(2, 1) += coef * (2.0 * (q1[0] * q2[1] + q1[2] * q2[3]));
}
}  // namespace

template <size_t MaxDeriv>
void UpdateRotationMatrices(
    gsl::not_null<std::array<Matrix, MaxDeriv + 1>*> rotation_matrices,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        fot_list,
    const std::string fot_name) {
  const auto& fot_ptr = fot_list.at(fot_name).get();

  if constexpr (MaxDeriv == 0) {
    std::array<DataVector, 1> quat = fot_ptr->func(t);
    Matrix R{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    AddBilinearTerm(R, quat[0], quat[0]);

    (*rotation_matrices)[0] = R;
  } else if constexpr (MaxDeriv == 1) {
    std::array<DataVector, 2> quat_and_deriv = fot_ptr->func_and_deriv(t);
    Matrix R{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    Matrix dtR{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    AddBilinearTerm(R, quat_and_deriv[0], quat_and_deriv[0]);

    AddBilinearTerm(dtR, quat_and_deriv[0], quat_and_deriv[1]);
    AddBilinearTerm(dtR, quat_and_deriv[1], quat_and_deriv[0]);

    (*rotation_matrices)[0] = R;
    (*rotation_matrices)[1] = dtR;
  } else if constexpr (MaxDeriv == 2) {
    std::array<DataVector, 3> quat_and_derivs = fot_ptr->func_and_2_derivs(t);
    Matrix R{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    Matrix dtR{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    Matrix dt2R{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    AddBilinearTerm(R, quat_and_derivs[0], quat_and_derivs[0]);

    AddBilinearTerm(dtR, quat_and_derivs[0], quat_and_derivs[1]);
    AddBilinearTerm(dtR, quat_and_derivs[1], quat_and_derivs[0]);

    AddBilinearTerm(dt2R, quat_and_derivs[0], quat_and_derivs[2]);
    AddBilinearTerm(dt2R, quat_and_derivs[1], quat_and_derivs[1], 2.0);
    AddBilinearTerm(dt2R, quat_and_derivs[2], quat_and_derivs[0]);

    (*rotation_matrices)[0] = R;
    (*rotation_matrices)[1] = dtR;
    (*rotation_matrices)[2] = dt2R;
  } else {
    ERROR("Cannot compute rotation matrix for MaxDeriv = " << MaxDeriv);
  }
}

// Generate instantiations for MaxDeriv = (0, 1, 2)
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                           \
  template void UpdateRotationMatrices<DIM(data)>(                     \
      gsl::not_null<std::array<Matrix, DIM(data) + 1>*>, const double, \
      const std::unordered_map<                                        \
          std::string,                                                 \
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&,  \
      const std::string);

GENERATE_INSTANTIATIONS(INSTANTIATE, (0, 1, 2))

#undef INSTANTIATE
#undef DIM
