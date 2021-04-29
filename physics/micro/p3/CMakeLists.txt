
set(MICRO_SRC fortran/micro_p3.F90
              fortran/micro_p3_types.F90
              fortran/micro_p3_utils.F90
              fortran/scream_abortutils.F90
              fortran/wv_sat_scream.F90)

add_library(micro ${MICRO_SRC})
target_compile_features(micro PUBLIC cxx_std_14)
target_include_directories(micro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(${MICRO_SRC} PROPERTIES COMPILE_FLAGS  "-ffree-line-length-none -fdefault-real-8 -fdefault-double-8")

