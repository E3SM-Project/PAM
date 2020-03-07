


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

# REALLY SHOULD PUT advection.h INTO JUST ADVECTION EXECUTABLE, ETC.
add_library(dycore STATIC src/models/advection.h src/grids/topology.h src/grids/geometry.h src/exchange.h src/variable_sets.h src/fields.h src/fileio.h src/timesteppers/RKSimple.h src/util.h)

set_target_properties(dycore PROPERTIES LINKER_LANGUAGE CXX)


#set_source_files_properties(src/fileio.h PROPERTIES COMPILE_FLAGS "-I${PNETCDF_INCLUDE}")


include_directories(src/driver)
include_directories(src/grids)
include_directories(src/models)
include_directories(src/operators)
include_directories(src/operators/differential)
include_directories(src/operators/reconstructions)
include_directories(src/timesteppers)
include_directories(src)

add_executable(advection driver/driver.cpp)
set_target_properties(advection PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(advection yakl dycore pnetcdf)

#"-L${PNETCDF_PATH}/lib -lpnetcdf"

#if ("${ARCH}" STREQUAL "CUDA")
#  enable_language(CUDA)
#  set_source_files_properties(driver.cpp PROPERTIES LANGUAGE CUDA)
#endif()

#-I${PNETCDF_INCLUDE}

set_source_files_properties(driver.cpp PROPERTIES COMPILE_FLAGS "${EXE_DEFS} ${YAKL_CXX_FLAGS}")

# -I${PNETCDF_PATH}/include
