# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_ParallelInterpolationRuntime")

set(LIBRARY_SOURCES
  )

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  # Boost::boost
  # ControlSystem
  # CoordinateMaps
  # DataStructures
  # Domain
  # DomainCreators
  # DomainStructure
  # Events
  # GeneralRelativitySolutions
  # GrSurfaces
  # H5
  # Interpolation
  # IoTestHelpers
  # Logging
  # MathFunctions
  # Observer
  # ParallelInterpolation
  # RelativisticEulerSolutions
  # Spectral
  # SphericalHarmonics
  ParallelInterpolationRuntime
  Utilities
  )
