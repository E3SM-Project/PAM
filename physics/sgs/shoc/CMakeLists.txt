
if (${PAM_SCREAM_USE_CXX})

  add_library(sgs INTERFACE)
  target_include_directories(sgs INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(sgs INTERFACE ${CMAKE_CURRENT_BINARY_DIR})

else()

  set(SGS_SRC fortran/shoc.F90)
  add_library(sgs ${SGS_SRC})
  target_include_directories(sgs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(sgs PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(sgs scream_common)
  set_source_files_properties(${SGS_SRC} PROPERTIES COMPILE_FLAGS  "-g -ffree-line-length-none -fdefault-real-8 -fdefault-double-8 ${YAKL_F90_FLAGS}")

endif()
