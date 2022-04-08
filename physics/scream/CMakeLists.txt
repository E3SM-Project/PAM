
set(SCREAM_F90_MODULES "${CMAKE_CURRENT_BINARY_DIR}/scream_f90_modules")
separate_arguments(SCREAM_F90_FLAGS UNIX_COMMAND "${YAKL_F90_FLAGS}")
set(SCREAM_Fortran_FLAGS ${SCREAM_F90_FLAGS} CACHE STRING "" FORCE)
set(SCREAM_LIB_ONLY TRUE)
option(EKAT_ENABLE_TESTS OFF FORCE)

add_subdirectory(../../externals/EKAT ./ekat                )

MESSAGE(STATUS "${CMAKE_CURRENT_SOURCE_DIR}/../../externals/ekat/cmake")
set (EKAT_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../externals/EKAT/cmake)
list(APPEND CMAKE_MODULE_PATH
     ${EKAT_CMAKE_PATH}
     ${EKAT_CMAKE_PATH}/mpi
     ${EKAT_CMAKE_PATH}/pkg_build
)

include(EkatMpiUtils)
DisableMpiCxxBindings()

if (Kokkos_ENABLE_CUDA)
  if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    if (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
      string (CONCAT msg
          "Kokkos_ENALBE_CUDA_RELOCATABLE_DEVICE_CODE=ON, and Kokkos_ENALBE_DEBUG_BOUNDS_CHECK=ON.\n"
          "   -> Disabling bounds checks, to prevent internal compiler errors.")
      message(WARNING "${msg}")
      set (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK OFF CACHE BOOL "" FORCE)
    else()
      string (CONCAT msg
          "Kokkos_ENALBE_CUDA_RELOCATABLE_DEVICE_CODE=ON.\n"
          "   -> Disabling bounds checks, to prevent internal compiler errors.")
      message(STATUS "${msg}")
      set (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK OFF CACHE BOOL "" FORCE)
    endif()
  endif()
  include(EkatBuildKokkos)
  # Note: we need Kokkos_SOURCE_DIR to be defined *before* calling EkatSetNvccWrapper.
  EkatSetKokkosSourceDir()

  include (EkatSetNvccWrapper)
  EkatSetNvccWrapper()
endif()

option (Kokkos_ENABLE_SERIAL "" ON)

###############################################################
## FOR KOKKOS DEBUG; REMOVE LATER TODO: 
###############################################################
# set(CMAKE_BUILD_TYPE_ci "debug" CACHE STRING "")
# set(Kokkos_ENABLE_DEBUG TRUE CACHE BOOL "")
# set(Kokkos_ENABLE_DEBUG_BOUNDS_CHECK TRUE CACHE BOOL "")
# set(Kokkos_ENABLE_DEBUG_DUALVIEW_MODIFY_CHECK TRUE CACHE BOOL "")
###############################################################
###############################################################
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/physics/shoc)
add_subdirectory(share                ./scream_share        )
add_subdirectory(physics/share        ./scream_physics_share)
add_subdirectory(physics/p3           ./scream_p3           )
add_subdirectory(physics/shoc         ./scream_shoc         )
add_dependencies(p3   physics_share)
add_dependencies(shoc physics_share)
target_include_directories(shoc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../pam_core)
target_include_directories(shoc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/physics/shoc)
target_include_directories(p3   PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../pam_core)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_physics_share)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_f90_modules  )

set(MY_OPT_FLAGS -O3)
target_compile_options(kokkoscore       PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)
target_compile_options(kokkoscontainers PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)
target_compile_options(kokkoscore       PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${MY_OPT_FLAGS}>)
target_compile_options(kokkoscontainers PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${MY_OPT_FLAGS}>)
target_compile_options(ekat             PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)
target_compile_options(physics_share    PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)
target_compile_options(p3               PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)
target_compile_options(shoc             PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)
target_compile_options(scream_share     PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${MY_OPT_FLAGS}>)

