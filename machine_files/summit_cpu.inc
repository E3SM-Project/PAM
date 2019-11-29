
KOKKOS_PATH = /ccs/home/imn/kokkos
KOKKOS_SRC_PATH = ${KOKKOS_PATH}
SRC = driver.cpp
EXE = cloudFV

default: main

CXX = mpic++
LINK = ${CXX}
# CFLAGS = -DARRAY_DEBUG -O1 -g -I${PNETCDF_PATH}/include
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

KOKKOS_DEVICES = "OpenMP"
KOKKOS_ARCH = "Power9"

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

include ${KOKKOS_PATH}/Makefile.kokkos

