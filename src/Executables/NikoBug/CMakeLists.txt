# Distributed under the MIT License.
# See LICENSE.txt for details.

set(EXECUTABLE NikoBug)

set(EXECUTABLE_SOURCES
  NikoBug.cpp
  ${EXECUTABLE}.decl.h
  )

add_charm_module(${EXECUTABLE})

add_spectre_executable(
  ${EXECUTABLE}
  EXCLUDE_FROM_ALL
  ${EXECUTABLE_SOURCES}
  )

add_dependencies(
  ${EXECUTABLE}
  # module_GlobalCache
  module_NikoBug
  )

target_link_libraries(
  ${EXECUTABLE}
  PRIVATE
  Charmxx::main
  )

