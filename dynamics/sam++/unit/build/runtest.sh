#!/bin/bash

printf "Rebuilding\n\n"

make -j8

printf "Running 2-D tests\n\n"

printf "Running Fortran code\n\n"
cd fortran2d
rm -f fortran_output_000001.nc
./fortran2d

printf "Running C++ code\n\n"
cd ../cpp2d
rm -f cpp_output_000001.nc
./cpp2d

printf "Comparing results\n\n"
cd ..
python nccmp.py fortran2d/fortran_output_000001.nc cpp2d/cpp_output_000001.nc



printf "Running 3-D tests\n\n"

printf "Running Fortran code\n\n"
cd fortran3d
rm -f fortran_output_000001.nc
./fortran3d

printf "Running C++ code\n\n"
cd ../cpp3d
rm -f cpp_output_000001.nc
./cpp3d

printf "Comparing results\n\n"
cd ..
python nccmp.py fortran3d/fortran_output_000001.nc cpp3d/cpp_output_000001.nc


