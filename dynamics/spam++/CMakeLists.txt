


cmake_minimum_required(VERSION 3.0)
project(spam++ CXX)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS}")

include_directories(${PNETCDF_INCLUDE})
link_directories(${PNETCDF_LIB})
link_libraries(pnetcdf)

add_subdirectory(externals/YAKL)
include_directories(externals/YAKL)
include_directories(externals/cub)
include_directories(externals/hipCUB)
include_directories(externals/rocPRIM)

include_directories(src)
include_directories(src/driver)
include_directories(src/grids)
include_directories(src/models)
include_directories(src/operators)
include_directories(src/operators/differential)
include_directories(src/operators/reconstructions)
include_directories(src/operators/cross)
include_directories(src/timesteppers)

if ("${MODEL}" STREQUAL "advection")
  include_directories(src/models/advection)
  add_executable(advection driver/driver.cpp)
  set_target_properties(advection PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(advection yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")
endif()

if ("${MODEL}" STREQUAL "dualadvection")
  include_directories(src/models/dualadvection)
  add_executable(dualadvection driver/driver.cpp)
  set_target_properties(dualadvection PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(dualadvection yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")
endif()

if ("${MODEL}" STREQUAL "swe")
  include_directories(src/models/swe)
  add_executable(swe driver/driver.cpp)
  set_target_properties(swe PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(swe yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")
endif()

if ("${MODEL}" STREQUAL "tswe")
  include_directories(src/models/tswe)
  add_executable(tswe driver/driver.cpp)
  set_target_properties(tswe PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(tswe yakl pnetcdf)
  set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")
endif()
