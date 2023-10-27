# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  RegisterDerivedWithCharm.cpp
  Sphere.cpp
)

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  RegisterDerivedWithCharm.hpp
  ShapeMapTransitionFunction.hpp
  Sphere.hpp
)
