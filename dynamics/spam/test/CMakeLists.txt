cmake_minimum_required(VERSION 3.0)
project(spam_tests)
enable_testing()

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

# YAKL
set(YAKL_HOME "${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/YAKL")
set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
add_subdirectory(${YAKL_HOME} ${YAKL_BIN})
include_directories(${YAKL_BIN})

# SPAM
set(PAMC_IO "none")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../ ./dycore)

# PAM core
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../pam_core pam_core)

include(${YAKL_HOME}/yakl_utils.cmake)

# Test directories
add_subdirectory(hodge_star)
