#!/bin/bash

./cmakeclean.sh

cmake    \
  -DPAM_SCREAM_CXX_FLAGS="${PAM_SCREAM_CXX_FLAGS}"                \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DYAKL_HOME="${YAKL_HOME}"                                      \
  -DSCREAM_HOME="${SCREAM_HOME}"                                  \
  -DSCREAM_DOUBLE_PRECISION=ON                                    \
  -DCMAKE_CUDA_HOST_COMPILER="`which mpic++`"                     \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  -DPAM_STANDALONE=ON                                          \
  -DKokkos_SOURCE_DIR=/home/mwarusz/kokkos \
  ../../../physics/scream_cxx_interfaces

