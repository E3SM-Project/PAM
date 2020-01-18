
add_executable(cpp2d ../dmdf.F90 ../cpp_driver.F90 ${CPP_SRC})
target_link_libraries(cpp2d yakl ${NCFLAGS})
set_property(TARGET cpp2d APPEND PROPERTY COMPILE_FLAGS ${DEFS2D} )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/kokkos/core/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../../externals/YAKL)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../kokkos)

