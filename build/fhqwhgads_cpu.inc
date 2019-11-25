
SRC = driver.cpp
EXE = cloudFV

default: main

CXX = mpic++
LINK = mpic++
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

