#!/bin/bash

./cmakeclean.sh

cmake      \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS}"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS}"     \
  -DYAKL_OPENMP45_FLAGS="${YAKL_OPENMP45_FLAGS}" \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS}"           \
  -DYAKL_C_FLAGS="${YAKL_C_FLAGS}"               \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"           \
  -DNC_INCLUDE="${NC_INCLUDE}"                   \
  -DNC_LIBS="${NC_LIBS}"                         \
  -DYAKL_ARCH="${YAKL_ARCH}"                     \
  ..
