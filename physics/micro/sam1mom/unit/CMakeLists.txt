cmake_minimum_required(VERSION 3.0)
project(sam1mom_unit)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/YAKL)
set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
add_subdirectory(${YAKL_HOME} ${YAKL_BIN})

set(F90_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/cloud.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/ice_fall.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/micro.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/precip_fall.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/precip_init.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/precip_proc.F90
                ${CMAKE_CURRENT_SOURCE_DIR}/../fortran/sat.F90 )

set(CXX_SOURCES driver.cpp)

set(DRIVER_SRC ${F90_SOURCES} ${CXX_SOURCES})

# Main driver
add_executable(driver ${DRIVER_SRC})
include_directories(${YAKL_BIN})
include_directories(${NC_INCLUDE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(driver)

set_source_files_properties(${F90_SOURCES} PROPERTIES COMPILE_FLAGS "${YAKL_F90_FLAGS}")

target_link_libraries(driver ${NC_LIBS})


