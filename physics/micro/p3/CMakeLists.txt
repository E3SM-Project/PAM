
set(MICRO_SRC micro_p3.F90
              micro_p3_interface.F90
              micro_p3_types.F90
              micro_p3_utils.F90
              scream_abortutils.F90
              wv_sat_scream.F90)

add_library(micro ${MICRO_SRC})
target_compile_features(micro PUBLIC cxx_std_14)
target_include_directories(micro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

