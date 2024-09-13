// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "ParallelAlgorithms/EventsAndTriggers/Completion.hpp"

namespace Events {
std::string Completion::name() const { return "Completion"; }

PUP::able::PUP_ID Completion::my_PUP_ID = 0;  // NOLINT
}  // namespace Events
