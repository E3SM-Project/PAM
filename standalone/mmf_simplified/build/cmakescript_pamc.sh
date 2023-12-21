#!/bin/bash

./cmakeclean.sh

# taken from https://gist.github.com/jrichardsz/c0047f58cb6765c7b6a7fb33c466ab3f#file-args_shell_parser-v1-0-0-md

PAM_MICRO="p3"
PAM_SGS="shoc"
PAMC_NDIMS="1"
PAMC_MODEL="extrudedmodel"
PAMC_HAMIL="man"
PAMC_THERMO="constkappavirpottemp"
PAMC_IO="serial"

for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done

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
  -DPAM_DYCORE="spam"                                             \
  -DPAM_MICRO=${PAM_MICRO}                                        \
  -DPAM_SGS=${PAM_SGS}                                            \
  -DPAM_RAD="none"                                                \
  -DPAMC_NDIMS=${PAMC_NDIMS}                                      \
  -DPAMC_MODEL=${PAMC_MODEL}                                      \
  -DPAMC_HAMIL=${PAMC_HAMIL}                                      \
  -DPAMC_THERMO=${PAMC_THERMO}                                    \
  -DPAMC_IO=${PAMC_IO}                                            \
  -DPAM_NLEV=${PAM_NLEV}                                          \
  -DSCREAM_CXX_LIBS_DIR=${SCREAM_CXX_LIBS_DIR}                    \
  -DPAM_SCREAM_USE_CXX=${PAM_SCREAM_USE_CXX}                      \
  ..
