
add_executable(cpp3d ../dmdf.F90 ../cpp_driver.F90 ${CPP_SRC})
target_link_libraries(cpp3d kokkos ${NCFLAGS})
set_property(TARGET cpp3d APPEND PROPERTY COMPILE_FLAGS ${DEFS3D} )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/kokkos/core/src)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../kokkos)

