
set(SGS_SRC fortran/shoc.F90)

add_library(sgs ${SGS_SRC})
target_compile_features(sgs PUBLIC cxx_std_14)
target_include_directories(sgs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(sgs PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(sgs scream_common physics_share)
set_source_files_properties(${SGS_SRC} PROPERTIES COMPILE_FLAGS  "-ffree-line-length-none -fdefault-real-8 -fdefault-double-8 ${YAKL_F90_FLAGS}")

