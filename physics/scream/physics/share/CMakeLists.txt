set(PHYSICS_SHARE_SRCS
          micro_p3_utils.F90
          physics_saturation.cpp
          physics_share.cpp
          physics_share_f2c.F90
          physics_test_data.cpp
          scream_abortutils.F90
          scream_types.F90
          wv_sat_scream.F90
)




# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND PHYSICS_SHARE_SRCS
  physics_saturation.cpp)
endif()

add_library(physics_share ${PHYSICS_SHARE_SRCS})
set_target_properties(physics_share PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_include_directories(physics_share PUBLIC ${SCREAM_F90_MODULES})
target_link_libraries(physics_share scream_share)
target_compile_options(physics_share PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

target_include_directories(physics_share PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(physics_share PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
