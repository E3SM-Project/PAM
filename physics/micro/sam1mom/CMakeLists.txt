
set(MICRO_SRC fortran/cloud.F90
              fortran/ice_fall.F90
              fortran/micro.F90
              fortran/precip_fall.F90
              fortran/precip_init.F90
              fortran/precip_proc.F90
              fortran/sat.F90 )

add_library(micro ${MICRO_SRC})
target_compile_features(micro PUBLIC cxx_std_14)
target_include_directories(micro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(${MICRO_SRC} PROPERTIES COMPILE_FLAGS  "-ffree-line-length-none -fdefault-real-8 -fdefault-double-8 ${YAKL_F90_FLAGS}")

