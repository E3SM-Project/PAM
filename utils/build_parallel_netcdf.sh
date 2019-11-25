#!/bin/bash

# This is a sample GNU parallel-netcdf build script
# Note that key environment variables are often not preserved in sudo
# Therefore, if placing in root folders, you should sudo this entire file

INSTALL_DIR=$HOME/parallel-netcdf-gnu

export MPICC=mpicc
export MPICXX=mpicxx
export MPIF77=mpif77
export MPIF90=mpif90
export CC=gcc
export CFLAGS='-fPIC -O2'
export FC=gfortran
export FCFLAGS='-fPIC -O2'
export F77=gfortran
export FFLAGS='-fPIC -O2'

make clean
make distclean
./configure --prefix=$INSTALL_DIR
make -j8
make install


