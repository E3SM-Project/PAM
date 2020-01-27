#!/bin/bash

./cmakeclean.sh

cmake \
  -DCMAKE_CXX_COMPILER=mpic++  \
  -DCMAKE_CXX_FLAGS="-O3"      \
  ..
