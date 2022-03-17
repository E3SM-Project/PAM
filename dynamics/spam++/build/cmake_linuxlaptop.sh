
#!/bin/bash

./build/cmake_clean.sh

cmake      \
  -DCMAKE_CXX_COMPILER=mpicxx                   \
  -DCXXFLAGS="-O3 -std=c++11"                   \
  -DPNETCDF_INCLUDE=/home/celdred/parallel-netcdf-gnu/include   \
  -DPNETCDF_LIB=/home/celdred/parallel-netcdf-gnu/lib   \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS} -DPAM_STANDALONE"     \
  -DYAKL_OPENMP45_FLAGS="${YAKL_OPENMP45_FLAGS} -DPAM_STANDALONE" \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_C_FLAGS="${YAKL_C_FLAGS}"                                \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"                            \
  -DNCFLAGS="${NCFLAGS}"                                          \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DMODEL=$1 -DGENERAL_COMPILE_CONST=$2 -DMODEL_COMPILE_CONST=$3
  ..
