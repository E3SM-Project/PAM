#!/bin/bash


############################################################################
## MAKE SURE WE HAVE WHAT WE NEED
############################################################################
function usage {
  printf "Usage: ./cmakescript.sh file.nc\n\n"
  printf "You must specify NCHOME and NFHOME environment variables to specify\n"
  printf "where the NetCDF libraries are located\n\n"
  printf "NetCDF binaries must include ncdump, nf-config, and nc-config\n\n"
  printf "You can also define FFLAGS to control optimizations and NCRMS \n"
  printf "to reduce the number of CRM samples and the runtime of the tests.\n\n"
}
if [[ "$1" == "" ]]; then
  printf "Error: missing NetCDF File parameter\n"
  usage
  exit -1
fi
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 0
fi
if [[ "$NCHOME" == "" ]]; then
  printf "Error: NCHOME environment variable not set\n"
  printf "set NCHOME with the path to the NetCDF C installation\n\n"
  usage
  exit -1
fi
if [[ "$NFHOME" == "" ]]; then
  printf "Error: NFHOME environment variable not set\n"
  printf "set NFHOME with the path to the NetCDF Fortran installation\n\n"
  usage
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

printf "Running with crm_nx=$NX, crm_ny=$NY, crm_nz=$NZ, crm_nx_rad=$NX_RAD, crm_ny_rad=$NY_RAD, crm_dx=$DX, crm_dt=$DT, yes3d=$YES3D, plev(nlev)=$PLEV, and infile=$INFILE\n\n"


############################################################################
## CLEAN UP THE PREVIOUS BUILD
############################################################################
rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake cpp CTestTestfile.cmake fortran Makefile fortran.exe cpp.exe


############################################################################
## GET THE NETCDF LINKING FLAGS
############################################################################
NCFLAGS="`$NFHOME/bin/nf-config --flibs` `$NCHOME/bin/nc-config --libs`"
printf "NetCDF Flags: $NCFLAGS\n\n"


############################################################################
## RUN THE CONFIGURE
############################################################################
FFLAGS="$FFLAGS -ffree-line-length-none -DCRM -DCRM_NX=$NX -DCRM_NY=$NY -DCRM_NZ=$NZ -DCRM_NX_RAD=$NX_RAD -DCRM_NY_RAD=$NY_RAD"
FFLAGS="$FFLAGS -DCRM_DT=$DT -DCRM_DX=$DX -DYES3DVAL=$YES3D -DPLEV=$PLEV -Dsam1mom -DINPUT_FILE=$INFILE"
FFLAGS="$FFLAGS -I$NCHOME/include -I$NFHOME/include  -O3"
if [[ "$NCRMS" != "" ]]; then
  FFLAGS="$FFLAGS -DNCRMS=$NCRMS"
fi
printf "FFLAGS: $FFLAGS\n\n"
cmake      \
    -DCMAKE_Fortran_FLAGS:STRING="$FFLAGS" \
    -DNCFLAGS:STRING="$NCFLAGS" \
    ..


