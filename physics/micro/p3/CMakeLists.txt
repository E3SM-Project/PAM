
set(MICRO_SRC fortran/micro_p3.F90)

add_library(micro ${MICRO_SRC})
target_include_directories(micro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(micro PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

message("USE_SCREAM is ${USE_SCREAM}")

target_link_libraries(micro scream_common)
set_source_files_properties(${MICRO_SRC} PROPERTIES COMPILE_FLAGS  "-ffree-line-length-none -fdefault-real-8 -fdefault-double-8 ${YAKL_F90_FLAGS}")

