# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_Ringdown")

set(LIBRARY_SOURCES
  Test_StrahlkorperCoefsInRingdownDistortedFrame.cpp
  )

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  Boost::boost
  CoordinateMaps
  DataStructures
  DomainCreators
  Domain
  DomainStructure
  ErrorHandling
  Options
  Ringdown
  Serialization
  Spectral
  SphericalHarmonics
  SphericalHarmonicsIO
  Utilities
  )
