#!/bin/bash -x

./cmakeclean.sh

if [[ "$1" != "" ]]; then
  . $1
fi

cmake      \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS} -DPAM_STANDALONE"     \
  -DYAKL_OPENMP45_FLAGS="${YAKL_OPENMP45_FLAGS} -DPAM_STANDALONE" \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_C_FLAGS="${YAKL_C_FLAGS}"                                \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"                            \
  -DPAM_LINK_FLAGS="${PAM_LINK_FLAGS}"                            \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DPAM_DYCORE="spam++"                                           \
  -DPAM_MICRO="p3"                                                \
  -DPAM_SGS="shoc"                                                \
  ${add_cmake_vars} \
  ..
