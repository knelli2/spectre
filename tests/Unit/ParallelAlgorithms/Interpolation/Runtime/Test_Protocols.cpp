// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "ParallelAlgorithms/Interpolation/Runtime/Points/SpecifiedPoints.hpp"
#include "ParallelAlgorithms/Interpolation/Runtime/Protocols/Points.hpp"
#include "Utilities/ProtocolHelpers.hpp"

static_assert(tt::assert_conforms_to_v<intrp2::points::SpecifiedPoints<3>,
                                       intrp2::protocols::Points>);
