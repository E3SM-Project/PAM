cmake_minimum_required(VERSION 3.0)
project(pam_standalone)

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

# Main driver
add_executable(driver driver.cpp)
include_directories(${YAKL_BIN})

#Choose dycore
if ("${PAM_DYCORE}" STREQUAL "AWFL")
  add_subdirectory(../dynamics/awfl ./dynamics)
endif()

#Choose microphysics
if ("${PAM_MICRO}" STREQUAL "kessler")
  add_subdirectory(../physics/micro/kessler ./micro)
endif()

target_link_libraries(driver dycore micro yakl ${NCFLAGS} -lyaml-cpp)

set_source_files_properties(${DRIVER_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${DRIVER_SRC} PROPERTIES LANGUAGE CUDA)
endif()


