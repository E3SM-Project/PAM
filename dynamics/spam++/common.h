#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include "Array.h"
#include "SArray.h"
#if defined(__USE_HIP__)
  #include "hip/hip_runtime.h"
#endif

using yakl::SArray;

typedef unsigned long ulong;
typedef unsigned int  uint;

// Define a real array type, presumed on the "device"
#if defined(__USE_CUDA__) || defined(__USE_HIP__)
  typedef yakl::Array<real,yakl::memDevice> realArr;
#else
  typedef yakl::Array<real,yakl::memHost> realArr;
#endif

// Define a real array type on the "host"
typedef yakl::Array<real,yakl::memHost> realArrHost;

// This will convert what would otherwise be double literals to the real type to avoid
// double precision arithmetic if real is single precision.
YAKL_INLINE real constexpr operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

// Specifying templated min and max functions
template <class T> YAKL_INLINE T mymin( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> YAKL_INLINE T mymax( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

// Declaring the precision for the model
typedef float real;

// Time scheme types
int constexpr TIME_TYPE_KGRK = 1;
int constexpr TIME_TYPE_ADER = 2;

// Reconstruction types

int constexpr RECONSTRUCTION_TYPE_FV   = 1;
int constexpr RECONSTRUCTION_TYPE_WENO = 2;








// COMPILE TIME CONSTANTS //

// Spatial order of accuracy for the model
int constexpr reconstruction_order_x = 1;
int constexpr reconstruction_order_y = 1;
int constexpr reconstruction_order_z = 1;
int constexpr differential_order_x = 2;
int constexpr differential_order_y = 2;
int constexpr differential_order_z = 2;

// Reconstruction type
int constexpr reconstruction_type_x = RECONSTRUCTION_TYPE_FV;
int constexpr reconstruction_type_y = RECONSTRUCTION_TYPE_FV;
int constexpr reconstruction_type_z = RECONSTRUCTION_TYPE_FV;

// Halo sizes
int constexpr maxhalosize_x = mymax(reconstruction_order_x,differential_order_x)/2; //IS THIS ALWAYS CORRECT?
int constexpr maxhalosize_z = mymax(reconstruction_order_y,differential_order_y)/2; //IS THIS ALWAYS CORRECT?
int constexpr maxhalosize_y = mymax(reconstruction_order_z,differential_order_z)/2; //IS THIS ALWAYS CORRECT?

// initial condition quadrature pts
int constexpr ic_quad_pts_x = 3;
int constexpr ic_quad_pts_y = 3;
int constexpr ic_quad_pts_z = 3;

// Time order of accuracy for the model
int constexpr time_order = 2;

// Time scheme
int constexpr time_type = TIME_TYPE_KGRK;

// Number of Dimensions
int constexpr ndims = 1;

#endif
