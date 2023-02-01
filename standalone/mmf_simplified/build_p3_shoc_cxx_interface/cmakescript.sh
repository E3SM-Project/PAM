#!/bin/bash

./cmakeclean.sh

cmake      \
  -DPAM_SCREAM_CXX_LINK_FLAGS="${PAM_SCREAM_CXX_LINK_FLAGS}"      \
  -DPAM_SCREAM_CXX_FLAGS="${PAM_SCREAM_CXX_FLAGS}"                \
  -DYAKL_HOME="${YAKL_HOME}"                                      \
  -DSCREAM_HOME="${SCREAM_HOME}"                                  \
  -DSCREAM_DOUBLE_PRECISION=ON                                    \
  -DCMAKE_CUDA_HOST_COMPILER="mpic++"                             \
  ../../../physics/scream_cxx_interfaces

