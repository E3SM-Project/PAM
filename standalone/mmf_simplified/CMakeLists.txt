cmake_minimum_required(VERSION 3.0)
project(pam_standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/YAKL)
set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
include_directories(${YAKL_BIN})

if (NOT ${PAM_SCREAM_USE_CXX})
  set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../externals/yaml-cpp ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../externals/eigen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../physics/scream_cxx_interfaces)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../pam_core pam_core)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../physics  physics )
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../dynamics dynamics)

add_executable(driver driver.cpp supercell_init.cpp supercell_init.F90)
add_executable(fortran_test fortran_test.F90)

target_compile_definitions(driver PUBLIC -DPAM_NLEV=${PAM_NLEV})

set_source_files_properties(supercell_init.F90 PROPERTY COMPILE_FLAGS "${YAKL_F90_FLAGS}")

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(driver)

if (PAM_SCREAM_USE_CXX)
  message(STATUS "*** building WITH C++ scream interface ***")
  target_compile_definitions(driver PUBLIC "-DP3_CXX -DSHOC_CXX")
  target_link_libraries(driver dynamics physics pam_core yaml-cpp "-L${CMAKE_CURRENT_BINARY_DIR}/../build_p3_shoc_cxx_interface -lpam_scream_cxx_interfaces -L${SCREAM_CXX_LIBS_DIR}/p3 -lp3 -L${SCREAM_CXX_LIBS_DIR}/shoc -lshoc -L${SCREAM_CXX_LIBS_DIR}/physics_share -lphysics_share -L${SCREAM_CXX_LIBS_DIR}/scream_share -lscream_share -L${SCREAM_CXX_LIBS_DIR}/externals/spdlog -lspdlog -L${SCREAM_CXX_LIBS_DIR} -lgptl -L${SCREAM_CXX_LIBS_DIR}/ekat/src/ekat -lekat -L${SCREAM_CXX_LIBS_DIR}/externals/kokkos/containers/src -lkokkoscontainers -L${SCREAM_CXX_LIBS_DIR}/externals/kokkos/core/src -lkokkoscore ${PAM_LINK_FLAGS} -ldl")
else()
  message(STATUS "*** building withOUT C++ scream interface ***")
  target_link_libraries(driver dynamics physics pam_core yaml-cpp ${PAM_LINK_FLAGS})
endif()

target_link_libraries(fortran_test dynamics physics pam_core ${PAM_LINK_FLAGS})

