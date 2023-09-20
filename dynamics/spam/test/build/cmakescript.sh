#!/bin/bash -x

./cmakeclean.sh

cmake      \
  -DCMAKE_CUDA_HOST_COMPILER=${CXX}              \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS}"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS}"     \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS}"           \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"           \
  -DPAM_LINK_FLAGS="${PAM_LINK_FLAGS}"           \
  -DYAKL_ARCH="${YAKL_ARCH}"                     \
  ..
