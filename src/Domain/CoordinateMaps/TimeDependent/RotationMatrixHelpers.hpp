// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines rotation matrices using quaternions.

#pragma once

#include <array>

#include "DataStructures/Matrix.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"

template <size_t MaxDeriv>
void UpdateRotationMatrices(
    gsl::not_null<std::array<Matrix, MaxDeriv + 1>*> rotation_matrices,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        fot_list,
    const std::string fot_name);
