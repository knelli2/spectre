# Distributed under the MIT License.
# See LICENSE.txt for details.

option(ASAN "Add AddressSanitizer compile flags" OFF)

# We handle the sanitizers using targets for the compile options, but modify
# CMAKE_EXE_LINKER_FLAGS for the linker flags because the sanitizers should
# only be linked into the final executable and CMake doesn't support linker
# interface flags before CMake 3.13.
add_library(Sanitizers IMPORTED INTERFACE)
add_library(Sanitizers::Address IMPORTED INTERFACE)
add_library(Sanitizers::Thread IMPORTED INTERFACE)
add_library(Sanitizers::UbInteger IMPORTED INTERFACE)
add_library(Sanitizers::UbUndefined IMPORTED INTERFACE)

if (ASAN)
  set_property(
    TARGET Sanitizers::Address
    APPEND PROPERTY
    INTERFACE_COMPILE_OPTIONS
    $<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer -fsanitize=address>
    )
  set(
    CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address"
    )
endif ()

option(THREADSAN "Add thread sanitizer behavior compile flags" OFF)
if (ASAN)
  set_property(
    TARGET Sanitizers::Thread
    APPEND PROPERTY
    INTERFACE_COMPILE_OPTIONS
    $<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer -fsanitize=thread>
    )
  set(
    CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=thread"
    )
  set(
    CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=thread"
    )
endif ()

option(UBSAN_UNDEFINED "Add UBSan undefined behavior compile flags" OFF)
if (UBSAN_UNDEFINED)
  set_property(
    TARGET Sanitizers::UbUndefined
    APPEND PROPERTY
    INTERFACE_COMPILE_OPTIONS
    $<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer -fsanitize=undefined>
    )
  set(
    CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=undefined"
    )
endif()

option(UBSAN_INTEGER "Add UBSan unsigned integer overflow compile flags" OFF)
if (UBSAN_INTEGER)
  set_property(
    TARGET Sanitizers::UbInteger
    APPEND PROPERTY
    INTERFACE_COMPILE_OPTIONS
    $<$<COMPILE_LANGUAGE:CXX>:-fno-omit-frame-pointer -fsanitize=integer>
    )
  set(
    CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer  -fsanitize=integer"
    )
  set(
    CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=integer"
    )
endif()

target_link_libraries(
  SpectreFlags
  INTERFACE
  Sanitizers::Address
  Sanitizers::Thread
  Sanitizers::UbInteger
  Sanitizers::UbUndefined
  )
