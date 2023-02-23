cmake_minimum_required(VERSION 3.8)
project(scream_cxx_interfaces CXX Fortran C)

set(MYSOURCES scream_cxx_interface_p3.cpp
              scream_cxx_interface_shoc.cpp
              scream_cxx_interface_finalize.cpp)

set(SCREAM_LIB_ONLY ON CACHE BOOL "" FORCE)
set(CMAKE_MODULE_PATH ${SCREAM_HOME}/components/eamxx/cmake;${SCREAM_HOME}/externals/ekat/cmake)
set(SCREAM_CIME_BUILD ON)
set(SCREAM_BASE_DIR ${SCREAM_HOME}/components/eamxx})
set(EKAT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ekat)
set(EKAT_DISABLE_TPL_WARNINGS ON)
set(EKAT_ENABLE_TESTS OFF)
set(Kokkos_ENABLE_SERIAL  ON CACHE BOOL "" FORCE)
# set(Kokkos_ENABLE_CUDA    ON CACHE BOOL "" FORCE)
# set(Kokkos_ENABLE_THREADS ON CACHE BOOL "" FORCE)
# set(Kokkos_ENABLE_HIP     ON CACHE BOOL "" FORCE)

set(SCREAM_POSSIBLY_NO_PACK_SIZE   FALSE       CACHE BOOL "" FORCE)
set(SCREAM_PACK_SIZE               1           CACHE BOOL "" FORCE)
set(SCREAM_SMALL_PACK_SIZE         1           CACHE BOOL "" FORCE)
set(SCREAM_NUM_VERTICAL_LEV        ${PAM_NLEV} CACHE BOOL "" FORCE)

include(EkatUtils)
EkatConfigFile(${SCREAM_HOME}/components/eamxx/src/scream_config.h.in
               ${SCREAM_HOME}/components/eamxx/src/scream_config.h
               F90_FILE ${SCREAM_HOME}/components/eamxx/src/scream_config.f)

include_directories(${SCREAM_HOME}/components/eamxx/src/physics/shoc )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/p3   )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/share)
include_directories(${SCREAM_HOME}/components/eamxx/src              )
include_directories(${YAKL_HOME}/external                            )

add_subdirectory(${SCREAM_HOME}/externals/ekat ${EKAT_BINARY_DIR})
add_library(pam_scream_cxx_interfaces ${MYSOURCES})
set_target_properties(pam_scream_cxx_interfaces PROPERTIES COMPILE_FLAGS ${PAM_SCREAM_CXX_FLAGS})
target_link_libraries(pam_scream_cxx_interfaces ekat)


