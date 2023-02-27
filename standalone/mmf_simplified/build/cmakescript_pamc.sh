#!/bin/bash

./cmakeclean.sh


export YAKL_CXX_FLAGS="-O3"

cmake      \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS} -DPAM_STANDALONE"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS} -DPAM_STANDALONE"     \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS} -DPAM_STANDALONE"           \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"                            \
  -DPAM_LINK_FLAGS="${PAM_LINK_FLAGS}"                            \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DCMAKE_CXX_COMPILER="mpic++"                                   \
  -DPAM_DYCORE="spam++"                                           \
  -DPAM_MICRO="p3"                                                \
  -DPAM_SGS="shoc"                                                \
  -DPAM_RAD="none"                                                \
  -DPAMC_MODEL="extrudedmodel"                                    \
  -DPAMC_HAMIL="man"                                              \
  -DPAMC_THERMO="constkappavirpottemp"                            \
  -DPAMC_IO="none"                                                \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  -DSCREAM_CXX_LIBS_DIR=${SCREAM_CXX_LIBS_DIR}                    \
  -DPAM_SCREAM_USE_CXX=${PAM_SCREAM_USE_CXX}                      \
  ..
