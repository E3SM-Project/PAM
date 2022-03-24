
add_library(dynamics INTERFACE)

#Choose dycore
if     ("${PAM_DYCORE}" STREQUAL "awfl")
  add_subdirectory(awfl ./dycore)
elseif ("${PAM_DYCORE}" STREQUAL "spam")
  message(FATAL_ERROR "ERROR: spam dycore doesn't exist yet")
  add_subdirectory(spam ./dycore)
else ()
  if ("${PAM_DYCORE}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_DYCORE=[awfl] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_DYCORE}\" is an invalid PAM_DYCORE CMake option. Valid options include \"awfl\".")
  endif()
endif()

target_link_libraries(dynamics INTERFACE dycore)

