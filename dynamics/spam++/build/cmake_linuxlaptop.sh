
#!/bin/bash

./build/cmake_clean.sh


cmake -DCMAKE_CXX_COMPILER=mpicxx                   \
      -DYAKL_CUB_HOME=`pwd`/../../../externals/cub   \
      -DCXXFLAGS="-O3 -std=c++11"                   \
      -DPNETCDF_INCLUDE=/home/celdred/parallel-netcdf-gnu/include   \
      -DPNETCDF_LIB=/home/celdred/parallel-netcdf-gnu/lib   \
      -DMODEL=$1 -DNTRACERS=$2 -DNTRACERS_FCT=$3 -DNTRACERS_Q=$4
