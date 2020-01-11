#!/bin/bash

############################################################################
## MAKE SURE WE HAVE WHAT WE NEED
############################################################################
if [[ "$1" == "" ]]; then
  echo "Error: missing NetCDF File parameter"
  echo "Usage: ./cmakescript.sh file.nc"
  exit -1
fi

if [[ "$NCHOME" == "" ]]; then
  echo "Error: NCHOME environment variable not set"
  echo "set NCHOME with the path to the NetCDF C installation"
  exit -1
fi

if [[ "$NFHOME" == "" ]]; then
  echo "Error: NFHOME environment variable not set"
  echo "set NFHOME with the path to the NetCDF Fortran installation"
  exit -1
fi


############################################################################
## GRAB DATA FROM THE NETCDF FILE
############################################################################
NX=`$NCHOME/bin/ncdump -h $1  | grep "crm_nx =" | awk '{print $3}'`
NY=`$NCHOME/bin/ncdump -h $1  | grep "crm_ny =" | awk '{print $3}'`
NZ=`$NCHOME/bin/ncdump -h $1  | grep "crm_nz =" | awk '{print $3}'`
NX_RAD=`$NCHOME/bin/ncdump -h $1  | grep "crm_nx_rad =" | awk '{print $3}'`
NY_RAD=`$NCHOME/bin/ncdump -h $1  | grep "crm_ny_rad =" | awk '{print $3}'`
DX=1000
DT=1
if [[ $NY -eq 1 ]]; then
  YES3D=0
else
  YES3D=1
fi
PLEV=`$NCHOME/bin/ncdump -h $1  | grep "nlev =" | awk '{print $3}'`
INFILE="\'$1\'"

echo "Running with crm_nx=$NX, crm_ny=$NY, crm_nz=$NZ, crm_nx_rad=$NX_RAD, crm_ny_rad=$NY_RAD, crm_dx=$DX, crm_dt=$DT, yes3d=$YES3D, plev(nlev)=$PLEV, and infile=$INFILE"


############################################################################
## CLEAN UP THE PREVIOUS BUILD
############################################################################
rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake cpp CTestTestfile.cmake fortran Makefile


############################################################################
## GET THE NETCDF LINKING FLAGS
############################################################################
NCFLAGS="`$NFHOME/bin/nf-config --flibs` `$NCHOME/bin/nc-config --libs`"
echo "NetCDF Flags: $NCFLAGS"


############################################################################
## RUN THE CONFIGURE
############################################################################
cmake      \
    -DCMAKE_Fortran_FLAGS:STRING=" -O3 -ffree-line-length-none -DCRM -DCRM_NX=$NX -DCRM_NY=$NY -DCRM_NZ=$NZ -DCRM_NX_RAD=$NX_RAD -DCRM_NY_RAD=$NY_RAD -DCRM_DT=$DT -DCRM_DX=$DX -DYES3DVAL=$YES3D -DPLEV=$PLEV -Dsam1mom -DINPUT_FILE=$INFILE -I$NCHOME/include -I$NFHOME/include " \
    -DNCFLAGS:STRING="$NCFLAGS" \
    ..


