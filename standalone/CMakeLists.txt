cmake_minimum_required(VERSION 3.0)
project(pam_standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../externals/YAKL)
set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
add_subdirectory(${YAKL_HOME} ${YAKL_BIN})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../externals/eigen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../coupler)

set(DRIVER_SRC driver.cpp)

# Main driver
add_executable(driver ${DRIVER_SRC})
include_directories(${YAKL_BIN})

#Choose dycore
if     ("${PAM_DYCORE}" STREQUAL "awfl")
  add_subdirectory(../dynamics/awfl ./dynamics)
elseif ("${PAM_DYCORE}" STREQUAL "spam")
  message(FATAL_ERROR "ERROR: spam dycore doesn't exist yet")
  # add_subdirectory(../dynamics/spam ./dynamics)
else ()
  if ("${PAM_DYCORE}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_DYCORE=[awfl] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_DYCORE}\" is an invalid PAM_DYCORE CMake option. Valid options include \"awfl\".")
  endif()
endif()

#Choose microphysics
if     ("${PAM_MICRO}" STREQUAL "kessler")
  add_subdirectory(../physics/micro/kessler ./micro)
elseif ("${PAM_MICRO}" STREQUAL "none")
  add_subdirectory(../physics/micro/none    ./micro)
elseif ("${PAM_MICRO}" STREQUAL "p3")
  add_subdirectory(../physics/micro/p3      ./micro)
else ()
  if ("${PAM_DYCORE}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_MICRO=[none|sat_adj|kessler|p3] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_DYCORE}\" is an invalid PAM_MICRO CMake option. Valid options include \"none\",\"sat_adj\",\"kessler\",\"p3\".")
  endif()
endif()

target_link_libraries(driver dycore micro yakl ${NCFLAGS} -lyaml-cpp)

if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${DRIVER_SRC} PROPERTIES LANGUAGE CUDA)
endif()
set_source_files_properties(${DRIVER_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS} -DPAM_STANDALONE")


