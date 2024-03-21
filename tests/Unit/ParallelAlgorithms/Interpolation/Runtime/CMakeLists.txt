# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_ParallelInterpolationRuntime")

set(LIBRARY_SOURCES
  Test_Protocols.cpp
  Test_Metafunctions.cpp
  )

add_subdirectory(Actions)
add_subdirectory(Callbacks)
add_subdirectory(Points)

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  # Boost::boost
  # ControlSystem
  # CoordinateMaps
  DataStructures
  # Domain
  # DomainCreators
  # DomainStructure
  # Events
  # GeneralRelativitySolutions
  GrSurfaces
  # H5
  # Interpolation
  # IoTestHelpers
  # Logging
  # MathFunctions
  # Observer
  # RelativisticEulerSolutions
  # Spectral
  # SphericalHarmonics
  ObserverHelpers
  Parallel
  ParallelInterpolationRuntime
  Utilities
  )
