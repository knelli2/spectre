# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY Ringdown)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  StrahlkorperCoefsInRingdownDistortedFrame.cpp
  WrapFillYlmData.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  StrahlkorperCoefsInRingdownDistortedFrame.hpp
  WrapFillYlmData.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  Boost::boost
  CoordinateMaps
  DataStructures
  DomainCreators
  Domain
  DomainStructure
  ErrorHandling
  Options
  ParallelInterpolation
  Serialization
  Spectral
  SphericalHarmonics
  SphericalHarmonicsIO
  Utilities
  )

add_subdirectory(Python)
