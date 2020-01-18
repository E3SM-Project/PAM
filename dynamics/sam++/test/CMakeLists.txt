cmake_minimum_required(VERSION 2.8)
project(standalone)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
enable_testing()

file(GLOB FORTRAN_SRC ../fortran/*.F90 ../fortran/SGS_TKE/*.F90 ../fortran/MICRO_SAM1MOM/*.F90 ../fortran/ADV_MPDATA/*.F90 ../fortran/*.c             )
file(GLOB CPP_SRC     ../cpp/*.F90     ../cpp/SGS_TKE/*.F90     ../cpp/MICRO_SAM1MOM/*.F90     ../cpp/ADV_MPDATA/*.F90     ../cpp/*.c     ../cpp/*.cpp)

add_subdirectory(../../../externals/YAKL   ./yakl  )

add_subdirectory(fortran2d)
add_subdirectory(fortran3d)
add_subdirectory(cpp2d)
add_subdirectory(cpp3d)

#add_test(NAME unit_abcoefs COMMAND unit_abcoefs)


