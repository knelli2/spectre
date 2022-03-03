# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  DarkEnergyFluid.cpp
  IdealFluid.cpp
  PolytropicFluid.cpp
  Spectral.cpp
  Enthalpy.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  DarkEnergyFluid.hpp
  EquationOfState.hpp
  IdealFluid.hpp
  PolytropicFluid.hpp
  Spectral.hpp
  Enthalpy.hpp
  )

add_subdirectory(Python)
