# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY ParallelInterpolationRuntime)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  )

add_dependencies(
  ${LIBRARY}
  module_GlobalCache
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  )

add_subdirectory(Callbacks)
add_subdirectory(Events)
add_subdirectory(Targets)
