
SRC = driver.cpp
EXE = cloudFV

default: main

CXX = nvcc -ccbin /usr/bin/mpic++ -expt-extended-lambda -x cu --std=c++11 --use_fast_math -expt-extended-lambda -arch=sm_35 -I/usr/local/cuda-10.1/include
LINK = nvcc -ccbin /usr/bin/mpic++
CXXFLAGS = -O3 -I${PNETCDF_PATH}/include
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

DEPFLAGS = -M
OBJ = $(notdir $(SRC:.cpp=.o))
LIB =

