#!/bin/bash

./cmakeclean.sh

unset GATOR_DISABLE

NETCDF_ROOT=/ccs/home/imn/software/netcdf_spock_gnu
YAML_ROOT=/ccs/home/imn/software/yaml_cpp_spock_gnu

export CC=hipcc
export CXX=hipcc
unset CXXFLAGS
export FFLAGS="-O3"

cmake -DYAKL_ARCH="HIP"                             \
      -DYAKL_HIP_FLAGS="-O3 -I${NETCDF_ROOT}/include -I${YAML_ROOT}/include -DYAKL_AUTO_PROFILE" \
      ..

