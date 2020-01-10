cmake_minimum_required(VERSION 2.8)
project(standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
enable_testing()

add_subdirectory(../fortran ./fortran)

add_executable(standalone driver.F90)
target_link_libraries(standalone fortran_code ${NCFLAGS})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/fortran)

#add_test(NAME unit_abcoefs COMMAND unit_abcoefs)


