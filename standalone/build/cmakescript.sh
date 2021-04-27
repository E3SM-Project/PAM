#!/bin/bash

./cmakeclean.sh

cmake      \
  -DCMAKE_CXX_FLAGS="${CXXFLAGS}"   \
  -DNCFLAGS="${NCFLAGS}"            \
  -DYAKL_HOME="${YAKL_HOME}"        \
  -DARCH="${ARCH}"                  \
  -DPAM_DYCORE="awfl"               \
  -DPAM_MICRO="p3"                  \
  ..


