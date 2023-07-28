
add_library(dynamics INTERFACE)

#Choose dycore
if     ("${PAM_DYCORE}" STREQUAL "awfl")
  add_subdirectory(awfl ./dycore)
  target_compile_definitions(dynamics INTERFACE PAM_DYCORE_AWFL)
elseif ("${PAM_DYCORE}" STREQUAL "spam++")
  add_subdirectory(spam ./dycore)
else ()
  if ("${PAM_DYCORE}" STREQUAL "")
    message(FATAL_ERROR "ERROR: You must specify -DPAM_DYCORE=[awfl,spam++] in your CMake configure line.")
  else ()
    message(FATAL_ERROR "ERROR: \"${PAM_DYCORE}\" is an invalid PAM_DYCORE CMake option. Valid options include \"awfl\",\"spam++\".")
  endif()
endif()

target_link_libraries(dynamics INTERFACE dycore)

