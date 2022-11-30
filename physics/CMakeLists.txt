

add_library(physics INTERFACE)

# Add scream common if it's needed
if ("${PAM_MICRO}" STREQUAL "p3" OR "${PAM_SGS}" STREQUAL "shoc")
  set(USE_SCREAM TRUE)
endif()

if (${USE_SCREAM})
  add_subdirectory(scream_common ./scream_common)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/scream_common)
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
target_link_libraries(physics INTERFACE yaml-cpp)

if (${USE_SCREAM})
  # cpp interface to kokkos version of scream
  target_include_directories(physics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cpp_interface)
endif()

