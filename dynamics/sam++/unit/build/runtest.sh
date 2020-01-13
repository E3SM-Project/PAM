#!/bin/bash

printf "Running 2-D tests\n\n"

./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc  || exit -1
make -j8                                                 || exit -1

rm -f cpp2d.nc
if [[ ! -f "fortran2d.nc" ]]; then
  ./fortran.exe                                          || exit -1
  mv fortran_output_000001.nc fortran2d.nc
fi
./cpp.exe                                                || exit -1
mv cpp_output_000001.nc cpp2d.nc

python nccmp.py cpp2d.nc fortran2d.nc                    || exit -1


echo "\n\n\nRunning 3-D tests\n\n"

./cmakescript.sh crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc   || exit -1
make -j8                                                 || exit -1

rm -f cpp3d.nc
if [[ ! -f "fortran3d.nc" ]]; then
  ./fortran.exe                                          || exit -1
  mv fortran_output_000001.nc fortran3d.nc
fi
./cpp.exe                                                || exit -1
mv cpp_output_000001.nc cpp3d.nc

python nccmp.py cpp2d.nc fortran2d.nc                    || exit -1

