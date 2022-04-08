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
  -DCMAKE_CUDA_HOST_COMPILER="mpic++"                             \
  -DPAM_DYCORE="awfl"                                             \
  -DPAM_MICRO="p3"                                                \
  -DPAM_SGS="shoc"                                                \
  -DKokkos_ENABLE_CUDA=ON                                         \
  -DKokkos_ARCH_KEPLER35=ON                                       \
  -DKokkos_ENABLE_CUDA_LAMBDA=ON                                  \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=ON                               \
  ..

