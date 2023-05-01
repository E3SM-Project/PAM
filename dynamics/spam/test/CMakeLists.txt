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
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/common")

# SPAM
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../ ./dycore)

# YAML-cpp
set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/yaml-cpp ${CMAKE_CURRENT_BINARY_DIR}/yaml-cpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/yaml-cpp/include)

# PAM core
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../../pam_core pam_core)

include(${YAKL_HOME}/yakl_utils.cmake)

macro(add_pamc_test testname model)
  add_executable(${testname} ${testname}.cpp)
  yakl_process_target(${testname})
  target_link_libraries(${testname} dycore yaml-cpp pam_core ${PAM_LINK_FLAGS})
  target_compile_definitions(${testname} PUBLIC ${model})
  add_test(NAME ${testname}_test COMMAND ./${testname})
endmacro(add_pamc_test)

# Test directories
add_subdirectory(hodge_star)
#add_subdirectory(fields)
#add_subdirectory(grid)
#add_subdirectory(hamiltonians)
#add_subdirectory(io)
#add_subdirectory(models)
add_subdirectory(operator_properties)
#add_subdirectory(parallel)
#add_subdirectory(timesteppers)
#
