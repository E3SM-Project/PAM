cmake_minimum_required(VERSION 2.8)
project(scampam)

add_executable(scampam driver/driver.cpp)

include_directories(${PNETCDF_INCLUDE})
link_directories(${PNETCDF_LIB})
link_libraries(pnetcdf)

add_subdirectory(externals/YAKL)
include_directories(externals/YAKL)
target_link_libraries(scampam yakl)

add_subdirectory(dynamics/awfl)
include_directories(dynamics/awfl)
target_link_libraries(scampam dycore)

