cmake_minimum_required(VERSION 3.8)
project(scream_cxx_interfaces CXX Fortran C)

set(MYSOURCES scream_cxx_interface_p3.cpp scream_cxx_interface_shoc.cpp)

set(SCREAM_LIB_ONLY ON CACHE BOOL "" FORCE)
set(CMAKE_MODULE_PATH ${SCREAM_HOME}/components/eamxx/cmake;${SCREAM_HOME}/externals/ekat/cmake)
set(SCREAM_CIME_BUILD ON)
set(SCREAM_BASE_DIR ${SCREAM_HOME}/components/eamxx})
set(EKAT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ekat)
set(EKAT_DISABLE_TPL_WARNINGS ON)
set(EKAT_ENABLE_TESTS OFF)
set(Kokkos_ENABLE_SERIAL  ON CACHE BOOL "" FORCE)
set(EKAT_CONFIGURE_FILE_F90_FILE ON)
# set(Kokkos_ENABLE_CUDA    ON CACHE BOOL "" FORCE)
# set(Kokkos_ENABLE_THREADS ON CACHE BOOL "" FORCE)
# set(Kokkos_ENABLE_HIP     ON CACHE BOOL "" FORCE)

set(SCREAM_POSSIBLY_NO_PACK_SIZE   FALSE CACHE BOOL "" FORCE)
set(SCREAM_PACK_SIZE               1     CACHE BOOL "" FORCE)
set(SCREAM_SMALL_PACK_SIZE         1     CACHE BOOL "" FORCE)
set(SCREAM_NUM_VERTICAL_LEV        50    CACHE BOOL "" FORCE)


macro (EkatConfigFile CONFIG_FILE_IN CONFIG_FILE_C EKAT_CONFIGURE_FILE_F90_FILE)
  message(STATUS "DEBUG: ${CONFIG_FILE_IN} ; ${CONFIG_FILE_C} ; ${EKAT_CONFIGURE_FILE_F90_FILE}")
  # Run the configure macro
  configure_file (${CONFIG_FILE_IN} ${CONFIG_FILE_C})
  # run sed to change '/*...*/' comments into '!/*...*/'
  execute_process(COMMAND cp ${CONFIG_FILE_C} ${EKAT_CONFIGURE_FILE_F90_FILE})
  execute_process(COMMAND sed -i "s;^/;!/;g" ${EKAT_CONFIGURE_FILE_F90_FILE})
  # # do the same for '//...' comments (turn them into '! ...'
  # execute_process(COMMAND sed 's;^//;!;g' ${EKAT_CONFIGURE_FILE_F90_FILE} > ${EKAT_CONFIGURE_FILE_F90_FILE}
  #                 WORKING_DIRECTORY ${EKAT_BINARY_DIR} 
  #                 ECHO_OUTPUT_VARIABLE ${TMP} )
endmacro (EkatConfigFile)


# include(EkatUtils)
EkatConfigFile(${SCREAM_HOME}/components/eamxx/src/scream_config.h.in
               ${SCREAM_HOME}/components/eamxx/src/scream_config.h
               ${SCREAM_HOME}/components/eamxx/src/scream_config.f)

include_directories(${SCREAM_HOME}/components/eamxx/src/physics/shoc )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/p3   )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/share)
include_directories(${SCREAM_HOME}/components/eamxx/src              )
include_directories(${SCREAM_HOME}/share/timing                      )
include_directories(${YAKL_HOME}/external                            )

add_definitions(-DSCREAM_CONFIG_IS_CMAKE)

add_subdirectory(${SCREAM_HOME}/externals/ekat                     ${EKAT_BINARY_DIR}                       )
add_subdirectory(scream_share                                      ${CMAKE_CURRENT_BINARY_DIR}/scream_share )
add_subdirectory(${SCREAM_HOME}/components/eamxx/src/physics/share ${CMAKE_CURRENT_BINARY_DIR}/physics_share)
add_subdirectory(${SCREAM_HOME}/components/eamxx/src/physics/p3    ${CMAKE_CURRENT_BINARY_DIR}/p3           )
add_subdirectory(${SCREAM_HOME}/components/eamxx/src/physics/shoc  ${CMAKE_CURRENT_BINARY_DIR}/shoc         )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(scream_share PROPERTIES COMPILE_FLAGS "-fpermissive -Wno-enum-compare")
endif()

