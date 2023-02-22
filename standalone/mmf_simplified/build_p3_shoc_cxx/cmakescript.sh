#!/bin/bash

./cmakeclean.sh

cmake      \
  -DPAM_SCREAM_CXX_FLAGS="${PAM_SCREAM_CXX_FLAGS}"                \
  -DYAKL_HOME="${YAKL_HOME}"                                      \
  -DSCREAM_HOME="${SCREAM_HOME}"                                  \
  -DSCREAM_DOUBLE_PRECISION=ON                                    \
  -DCMAKE_CUDA_HOST_COMPILER="mpic++"                             \
  -DSCREAM_Fortran_FLAGS="${SCREAM_Fortran_FLAGS}"                \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  ../../../physics/scream_cxx_p3_shoc

