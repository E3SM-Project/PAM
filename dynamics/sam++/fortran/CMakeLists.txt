
file(GLOB SRC *.F90 MICRO_SAM1MOM/*.F90 SGS_TKE/*.F90 ADV_MPDATA/*.F90)

add_library(fortran_code STATIC ${SRC})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

