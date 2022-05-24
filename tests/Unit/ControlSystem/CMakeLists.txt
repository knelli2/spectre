# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_ControlSystem")

set(LIBRARY_SOURCES
  )

add_subdirectory(Actions)
add_subdirectory(ControlErrors)
add_subdirectory(Protocols)
add_subdirectory(Systems)
add_subdirectory(Tags)

add_test_library(
  ${LIBRARY}
  "ControlSystem"
  "${LIBRARY_SOURCES}"
  "Boost::boost;ControlSystem;Domain;FunctionsOfTime;ApparentHorizons"
  )
