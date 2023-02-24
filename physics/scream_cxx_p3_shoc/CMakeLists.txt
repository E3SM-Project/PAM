cmake_minimum_required(VERSION 3.8)
project(scream_cxx_interfaces CXX Fortran C)

set(MYSOURCES scream_cxx_interface_p3.cpp
              scream_cxx_interface_shoc.cpp)


set(SCREAM_LIB_ONLY ON CACHE BOOL "" FORCE)
set(CMAKE_MODULE_PATH ${SCREAM_HOME}/components/eamxx/cmake;${SCREAM_HOME}/externals/ekat/cmake)
set(SCREAM_CIME_BUILD ON)
set(SCREAM_BASE_DIR ${SCREAM_HOME}/components/eamxx})
set(EKAT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/ekat)
option(EKAT_DISABLE_TPL_WARNINGS "Whether we should suppress warnings when compiling TPLs." ON)
option(EKAT_ENABLE_TESTS "Whether tests should be built." OFF)
set(EKAT_CONFIGURE_FILE_F90_FILE ON)
set(Kokkos_ENABLE_SERIAL  ON CACHE BOOL "" FORCE)
if ("${YAKL_ARCH}" STREQUAL "CUDA")
  set(EAMXX_ENABLE_GPU             ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_CUDA_LAMBDA    ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_CUDA           ON CACHE BOOL "" FORCE)
elseif ("${YAKL_ARCH}" STREQUAL "HIP")
  set(EAMXX_ENABLE_GPU                     ON    CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_HIP                    ON    CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_DEPRECATED_CODE        FALSE CACHE BOOL "")
  set(Kokkos_ENABLE_EXPLICIT_INSTANTIATION FALSE CACHE BOOL "")
  set($ENV{Kokkos_ARCH_HIP}                ON    CACHE BOOL "" FORCE)
  set(CMAKE_CXX_STANDARD 14)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang")
  set(CMAKE_Fortran_FLAGS "-Mpreprocess")
endif()

add_definitions(-DSCREAM_CIME_BUILD -DHAVE_MPI -DSCREAM_CONFIG_IS_CMAKE)

set(SCREAM_POSSIBLY_NO_PACK_SIZE   FALSE       CACHE BOOL "" FORCE)
set(SCREAM_PACK_SIZE               1           CACHE BOOL "" FORCE)
set(SCREAM_SMALL_PACK_SIZE         1           CACHE BOOL "" FORCE)
set(SCREAM_NUM_VERTICAL_LEV        ${PAM_NLEV} CACHE BOOL "" FORCE)
set(SCREAM_DATA_DIR                ${CMAKE_CURRENT_SOURCE_DIR}/../micro/p3)

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

EkatConfigFile(${SCREAM_HOME}/components/eamxx/src/scream_config.h.in
               ${SCREAM_HOME}/components/eamxx/src/scream_config.h
               ${SCREAM_HOME}/components/eamxx/src/scream_config.f)

include_directories(${SCREAM_HOME}/components/eamxx/src/physics/shoc     )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/shoc/impl)
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/p3       )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/p3/impl  )
include_directories(${SCREAM_HOME}/components/eamxx/src/physics/share    )
include_directories(${SCREAM_HOME}/components/eamxx/src                  )
include_directories(${SCREAM_HOME}/components/eam/src/physics/cam        )
include_directories(${SCREAM_HOME}/share/timing                          )
include_directories(${YAKL_HOME}/external                                )

add_subdirectory(${SCREAM_HOME}/externals/ekat                     ${EKAT_BINARY_DIR}                       )
add_subdirectory(scream_share                                      ${CMAKE_CURRENT_BINARY_DIR}/scream_share )
add_subdirectory(${SCREAM_HOME}/components/eamxx/src/physics/share ${CMAKE_CURRENT_BINARY_DIR}/physics_share)
add_subdirectory(${SCREAM_HOME}/components/eamxx/src/physics/p3    ${CMAKE_CURRENT_BINARY_DIR}/p3           )
add_subdirectory(${SCREAM_HOME}/components/eamxx/src/physics/shoc  ${CMAKE_CURRENT_BINARY_DIR}/shoc         )
add_library(gptl  ${SCREAM_HOME}/externals/scorpio/src/gptl/GPTLget_memusage.c
                  ${SCREAM_HOME}/externals/scorpio/src/gptl/GPTLprint_memusage.c
                  ${SCREAM_HOME}/externals/scorpio/src/gptl/GPTLutil.c
                  ${SCREAM_HOME}/externals/scorpio/src/gptl/f_wrappers.c
                  ${SCREAM_HOME}/externals/scorpio/src/gptl/gptl.c
                  ${SCREAM_HOME}/externals/scorpio/src/gptl/gptl_papi.c)


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(scream_share PROPERTIES COMPILE_FLAGS "-fpermissive -Wno-enum-compare")
endif()

