cmake_minimum_required(VERSION 3.0)
project(yakl_unit)
enable_testing()

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

# YAKL library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/YAKL ./yakl)

# Unit test directories
add_subdirectory(recon_regular)
add_subdirectory(recon_irregular)
add_subdirectory(anelastic_direct)

