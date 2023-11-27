#!/bin/bash

./cmakeclean.sh

if [[ "${YAKL_ARCH}" == "HIP" ]]; then
  sed -i 's/+bfb_erf/+ bfb_erf/' ${SCREAM_HOME}/components/eam/src/physics/cam/shoc.F90
fi

cmake      \
  -DYAKL_ARCH="${YAKL_ARCH}"                                      \
  -DYAKL_HOME="${YAKL_HOME}"                                      \
  -DSCREAM_HOME="${SCREAM_HOME}"                                  \
  -DSCREAM_DOUBLE_PRECISION=ON                                    \
  -DCMAKE_CUDA_HOST_COMPILER="`which mpic++`"                     \
  -DSCREAM_Fortran_FLAGS="${SCREAM_Fortran_FLAGS}"                \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  -DKokkos_SOURCE_DIR=/home/mwarusz/kokkos \
  -DPAM_STANDALONE=ON \
  ../../../physics/scream_cxx_p3_shoc

