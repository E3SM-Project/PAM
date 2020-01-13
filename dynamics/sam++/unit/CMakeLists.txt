cmake_minimum_required(VERSION 2.8)
project(standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
enable_testing()

add_subdirectory(../fortran ./fortran)
add_subdirectory(../cpp     ./cpp    )

add_executable(fortran.exe dmdf.F90 fortran_driver.F90)
target_link_libraries(fortran.exe fortran_code ${NCFLAGS})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/fortran)

add_executable(cpp.exe dmdf.F90 cpp_driver.F90)
target_link_libraries(cpp.exe cpp_code ${NCFLAGS})
include_directories(${CMAKE_CURRENT_BINARY_DIR}/cpp)

#add_test(NAME unit_abcoefs COMMAND unit_abcoefs)


