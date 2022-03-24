

add_library(physics INTERFACE)

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
  include_directories(scream)
  include_directories(scream/physics/shoc)
  set(SCREAM_F90_MODULES "${CMAKE_CURRENT_BINARY_DIR}/scream_f90_modules")
  separate_arguments(SCREAM_F90_FLAGS UNIX_COMMAND "${YAKL_F90_FLAGS}")
  set(SCREAM_Fortran_FLAGS ${SCREAM_F90_FLAGS} CACHE STRING "" FORCE)
  set(SCREAM_LIB_ONLY TRUE)
  option(EKAT_ENABLE_TESTS OFF FORCE)
  add_subdirectory(../externals/EKAT     ./ekat                )
  add_subdirectory(scream/share          ./scream_share        )
  add_subdirectory(scream/physics/share  ./scream_physics_share)
  add_subdirectory(scream_common         ./scream_common       )
  add_subdirectory(scream/physics/p3     ./scream_p3           )
  add_dependencies(scream_common         physics_share         )
  add_dependencies(p3                    physics_share         )
  if ("${PAM_SGS}" STREQUAL "shoc")
    add_subdirectory(scream/physics/shoc  ./scream_shoc)
    add_dependencies(shoc                 physics_share)
    target_include_directories(shoc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../pam_core)
    target_include_directories(shoc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/scream/physics/shoc)
  endif()
  target_include_directories(p3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../pam_core)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_common       )
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_physics_share)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_f90_modules  )
endif()

#Choose microphysics
if     ("${PAM_MICRO}" STREQUAL "kessler")
  add_subdirectory(micro/kessler ./micro)
elseif ("${PAM_MICRO}" STREQUAL "none")
  add_subdirectory(micro/none    ./micro)
elseif ("${PAM_MICRO}" STREQUAL "p3")
  add_subdirectory(micro/p3      ./micro)
elseif ("${PAM_MICRO}" STREQUAL "sam1mom")
  add_subdirectory(micro/sam1mom ./micro)
else ()
  if ("${PAM_MICRO}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_MICRO=[none|kessler|p3|sam1mom] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_MICRO}\" is an invalid PAM_MICRO CMake option. Valid options include \"none\",\"sam1mom\",\"kessler\",\"p3\".")
  endif()
endif()

#Choose SGS
if     ("${PAM_SGS}" STREQUAL "shoc")
  add_subdirectory(sgs/shoc ./sgs)
elseif ("${PAM_SGS}" STREQUAL "none")
  add_subdirectory(sgs/none ./sgs)
endif()

target_link_libraries(physics INTERFACE micro sgs)
if (${USE_SCREAM})
  target_link_libraries(physics INTERFACE scream_common p3)
  if ("${PAM_SGS}" STREQUAL "shoc")
    target_link_libraries(physics INTERFACE shoc)
  endif()
else()
  target_link_libraries(physics INTERFACE yaml-cpp)
endif()

