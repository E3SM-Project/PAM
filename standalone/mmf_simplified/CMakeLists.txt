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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../externals/eigen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../coupler)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../mmf_interface)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../mmf_routines)

set(DRIVER_SRC driver.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/../../mmf_interface/mmf_interface_globals.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/../../mmf_interface/mmf_interface_extern_c.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/../../mmf_interface/mmf_fortran_interface.F90
               ${CMAKE_CURRENT_SOURCE_DIR}/../../coupler/pam_coupler_globals.cpp)

# Main driver
add_executable(driver ${DRIVER_SRC})
include_directories(${YAKL_BIN})

#Choose dycore
if     ("${PAM_DYCORE}" STREQUAL "awfl")
  add_subdirectory(../../dynamics/awfl ./dynamics)
elseif ("${PAM_DYCORE}" STREQUAL "spam")
  message(FATAL_ERROR "ERROR: spam dycore doesn't exist yet")
  # add_subdirectory(../../dynamics/spam ./dynamics)
else ()
  if ("${PAM_DYCORE}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_DYCORE=[awfl] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_DYCORE}\" is an invalid PAM_DYCORE CMake option. Valid options include \"awfl\".")
  endif()
endif()

# Add scream common if it's needed
if ("${PAM_MICRO}" STREQUAL "p3" OR "${PAM_SGS}" STREQUAL "shoc")
  set(USE_SCREAM TRUE)
endif()

if (${USE_SCREAM})
  ###############################################################
  ## FOR KOKKOS DEBUG; REMOVE LATER TODO: 
  ###############################################################
  # set(CMAKE_BUILD_TYPE_ci "debug" CACHE STRING "")
  # set(Kokkos_ENABLE_DEBUG TRUE CACHE BOOL "")
  # set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK TRUE CACHE BOOL "")
  # set(Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK TRUE CACHE BOOL "")
  ###############################################################
  ###############################################################
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../physics/scream)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../physics/scream/physics/shoc)
  set(SCREAM_F90_MODULES "${CMAKE_CURRENT_BINARY_DIR}/scream_f90_modules")
  separate_arguments(SCREAM_F90_FLAGS UNIX_COMMAND "${YAKL_F90_FLAGS}")
  set(SCREAM_Fortran_FLAGS ${SCREAM_F90_FLAGS} CACHE STRING "" FORCE)
  set(SCREAM_LIB_ONLY TRUE)
  option(EKAT_ENABLE_TESTS OFF FORCE)
  add_subdirectory(../../externals/EKAT                ./ekat                )
  add_subdirectory(../../physics/scream/share          ./scream_share        )
  add_subdirectory(../../physics/scream/physics/share  ./scream_physics_share)
  add_subdirectory(../../physics/scream_common         ./scream_common       )
  add_dependencies(scream_common physics_share)
  add_subdirectory(../../physics/scream/physics/p3     ./scream_p3           )
  add_dependencies(p3 physics_share)
  if ("${PAM_SGS}" STREQUAL "shoc")
    add_subdirectory(../../physics/scream/physics/shoc   ./scream_shoc         )
    add_dependencies(shoc physics_share)
  endif()
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_common       )
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_physics_share)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_f90_modules  )
endif()

#Choose microphysics
if     ("${PAM_MICRO}" STREQUAL "kessler")
  add_subdirectory(../../physics/micro/kessler ./micro)
elseif ("${PAM_MICRO}" STREQUAL "none")
  add_subdirectory(../../physics/micro/none    ./micro)
elseif ("${PAM_MICRO}" STREQUAL "p3")
  add_subdirectory(../../physics/micro/p3      ./micro)
elseif ("${PAM_MICRO}" STREQUAL "wsm6")
  add_subdirectory(../../physics/micro/wsm6    ./micro)
elseif ("${PAM_MICRO}" STREQUAL "sam1mom")
  add_subdirectory(../../physics/micro/sam1mom ./micro)
else ()
  if ("${PAM_MICRO}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_MICRO=[none|kessler|p3|sam1mom] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_DYCORE}\" is an invalid PAM_MICRO CMake option. Valid options include \"none\",\"sam1mom\",\"kessler\",\"p3\".")
  endif()
endif()

#Choose SGS
if     ("${PAM_SGS}" STREQUAL "shoc")
  add_subdirectory(../../physics/sgs/shoc ./sgs)
elseif ("${PAM_SGS}" STREQUAL "none")
  add_subdirectory(../../physics/sgs/none ./sgs)
endif()

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(driver)

target_link_libraries(driver dycore micro sgs)
if (${USE_SCREAM})
  target_link_libraries(driver scream_common p3)
  if ("${PAM_SGS}" STREQUAL "shoc")
    target_link_libraries(driver shoc)
  endif()
endif()
target_link_libraries(driver ${NCFLAGS})


