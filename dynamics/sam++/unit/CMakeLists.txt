cmake_minimum_required(VERSION 2.8)
project (sam_unit CXX)

enable_language(Fortran)
enable_testing()

add_subdirectory(../fortran ./fortran)
add_subdirectory(../cpp ./cpp)

add_executable(unit_abcoefs unit_abcoefs.cpp)
target_link_libraries(unit_abcoefs fortran_code cpp_code)
include_directories(../cpp ../fortran)
add_test(NAME unit_abcoefs COMMAND unit_abcoefs)


