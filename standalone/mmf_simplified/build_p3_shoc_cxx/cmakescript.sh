#!/bin/bash

./cmakeclean.sh

# sed -i  's/GPTLpr_summary_file (comm.mpi_comm(),fname.c_str());/GPTLpr_summary_file ((int)comm.mpi_comm(),fname.c_str());/' ${SCREAM_HOME}/components/eamxx/src/share/util/scream_timing.cpp

cmake      \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DYAKL_HOME="${YAKL_HOME}"                                      \
  -DSCREAM_HOME="${SCREAM_HOME}"                                  \
  -DSCREAM_DOUBLE_PRECISION=ON                                    \
  -DCMAKE_CUDA_HOST_COMPILER="`which mpic++`"                     \
  -DSCREAM_Fortran_FLAGS="${SCREAM_Fortran_FLAGS}"                \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  ../../../physics/scream_cxx_p3_shoc

