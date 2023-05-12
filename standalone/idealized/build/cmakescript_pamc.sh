#!/bin/bash -x

./cmakeclean.sh

cmake      \
  -DCMAKE_CXX_COMPILER=${CXX}                                     \
  -DCMAKE_CUDA_HOST_COMPILER=${CXX}                               \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS} -DPAM_STANDALONE"     \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"                            \
  -DPAM_LINK_FLAGS="${PAM_LINK_FLAGS}"                            \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DPAM_DYCORE="spam++"                                           \
  -DPAM_MICRO="none"                                              \
  -DPAM_SGS="none"                                                \
  -DPAM_RAD="none"                                                \
  -DPAMC_MODEL="extrudedmodel"                                    \
  -DPAMC_HAMIL="an"                                              \
  -DPAMC_THERMO="idealgaspottemp"                            \
  -DPAMC_IO="serial"                                              \
  ..
