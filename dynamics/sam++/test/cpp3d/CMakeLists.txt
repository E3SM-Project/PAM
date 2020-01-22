
add_executable(cpp3d ../dmdf.F90 ../cpp_driver.F90 ${CPP_SRC})
target_link_libraries(cpp3d yakl ${NCFLAGS})
set_property(TARGET cpp3d APPEND PROPERTY COMPILE_FLAGS ${DEFS3D} )

if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${CUDA_SRC} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS " -D__AUTO_FENCE__-D__USE_CUDA__ -D__MANAGED__ -std=c++11 --expt-extended-lambda")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/YAKL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/cub)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/hipCUB)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/rocPRIM)

