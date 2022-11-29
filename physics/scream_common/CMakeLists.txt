
set(SCREAM_SRC scream_types.F90
               micro_p3_utils.F90
               scream_abortutils.F90
               wv_sat_scream.F90)

add_library(scream_common ${SCREAM_SRC})
target_include_directories(scream_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(scream_common PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(scream_common)
set_source_files_properties(${SCREAM_SRC} PROPERTIES COMPILE_FLAGS  "-ffree-line-length-none -fdefault-real-8 -fdefault-double-8 ${YAKL_F90_FLAGS}")

