
file(GLOB SRC *.F90 SGS_TKE/*.F90 MICRO_SAM1MOM/*.F90 ADV_MPDATA/*.F90 *.c)

add_library(fortran_code STATIC ${SRC})

