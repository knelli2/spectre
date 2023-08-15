// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"

namespace domain::FunctionsOfTime {
/*!
 * \brief Put time bounds for all functions of time in a nicely formatted string
 */
std::string ouput_time_bounds(
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);
}  // namespace domain::FunctionsOfTime
