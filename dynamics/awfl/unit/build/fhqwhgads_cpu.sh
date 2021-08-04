#!/bin/bash

./cmakeclean.sh

unset GATOR_DISABLE

export CC=gcc
export CXX=g++
export CXXFLAGS="-O3 -DYAKL_AUTO_PROFILE"
export FFLAGS="-O3"

cmake ..

