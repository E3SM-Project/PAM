
file(GLOB SRC *.F90 SGS_TKE/*.F90 MICRO_SAM1MOM/*.F90 ADV_MPDATA/*.F90 *.c *.cpp)

add_library(cpp_code STATIC ${SRC})

