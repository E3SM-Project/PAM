#!/bin/bash

./cmakeclean.sh

cmake      \
  -DCMAKE_CXX_COMPILER=${CXX}                                     \
  -DCMAKE_CUDA_HOST_COMPILER=${CXX}                               \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS} -DPAM_STANDALONE -DPAMA_DYCORE"         \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS} -DPAM_STANDALONE -DPAMA_DYCORE"           \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS} -DPAM_STANDALONE -DPAMA_DYCORE"         \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS} -DPAM_STANDALONE -DPAMA_DYCORE"     \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS} -DPAM_STANDALONE -DPAMA_DYCORE"           \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"                            \
  -DPAM_LINK_FLAGS="${PAM_LINK_FLAGS}"                            \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DPAM_DYCORE="awfl"                                             \
  -DPAM_MICRO="p3"                                                \
  -DPAM_RAD="none"                                                \
  -DPAM_SGS="none"                                                \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  -DSCREAM_CXX_LIBS_DIR=${SCREAM_CXX_LIBS_DIR}                    \
  -DPAM_SCREAM_USE_CXX=${PAM_SCREAM_USE_CXX}                      \
  ..

