
add_executable(cpp2d ../dmdf.F90 ../cpp_driver.F90 ${CPP_SRC})
target_link_libraries(cpp2d yakl ${NCFLAGS})
set_property(TARGET cpp2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )
set_source_files_properties(${CUDA_SRC} PROPERTIES LANGUAGE CUDA)
set_source_files_properties(${CUDA_SRC} PROPERTIES COMPILE_FLAGS "-D__USE_CUDA__ -D__MANAGED__ -std=c++11 --expt-extended-lambda")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/YAKL)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/cub)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/hipCUB)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/rocPRIM)

