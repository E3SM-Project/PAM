cmake_minimum_required(VERSION 3.0)
project(driver)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

set(DRIVER_SRC driver.cpp)

set(YAKL_HOME ${CMAKE_CURRENT_BINARY_DIR}/../../../externals/YAKL)
add_subdirectory(${YAKL_HOME} ./yakl)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../externals/eigen)

# Main driver
add_executable(driver ${DRIVER_SRC})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/yakl)
include_directories(${YAKL_HOME})
target_link_libraries(driver yakl ${NCFLAGS} -lyaml-cpp)

set_source_files_properties(${DRIVER_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${DRIVER_SRC} PROPERTIES LANGUAGE CUDA)
endif()

target_compile_features(driver PUBLIC cxx_std_14)

