#!/bin/bash

./cmakeclean.sh

cmake      \
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
  -DPAM_DYCORE="spam++"                                             \
  -DPAM_MICRO="none"                                           \
  -DPAM_SGS="none"                                                \
  -DPAMC_MODEL="extrudedmodel"                                                \
  -DPAMC_HAMIL="mce_rho"                                                \
  -DPAMC_THERMO="constkappavirpottemp"                                                \
  -DPAMC_IO="serial"                                                \
  ..
