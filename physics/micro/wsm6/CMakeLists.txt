
set(MICRO_SRC module_model_constants.F90
              module_mp_radar.F90
              module_mp_wsm6.F90
              module_wrf_error.F90 )

add_library(micro ${MICRO_SRC})
target_compile_features(micro PUBLIC cxx_std_14)
target_include_directories(micro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_source_files_properties(${MICRO_SRC} PROPERTIES COMPILE_FLAGS "-fdefault-real-8 -fdefault-double-8")

