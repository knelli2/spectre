# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY Ringdown)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  StrahlkorperCoefsInRingdownDistortedFrame.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  StrahlkorperCoefsInRingdownDistortedFrame.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  Boost::boost
  DataStructures
  DomainCreators
  Domain
  DomainStructure
  ErrorHandling
  Options
  Serialization
  Spectral
  SphericalHarmonics
  SphericalHarmonicsIO
  Utilities
  )

add_subdirectory(Python)
