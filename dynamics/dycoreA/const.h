
#ifndef _CONST_H_
#define _CONST_H_

#include <cmath>
#include "Array.h"
#include "SArray.h"

using yakl::SArray;

typedef float         real;
typedef unsigned long ulong;
typedef unsigned int  uint;

#if defined(__USE_CUDA__) || defined(__USE_HIP__)
typedef yakl::Array<real,yakl::memDevice> realArr;
#else
typedef yakl::Array<real,yakl::memHost> realArr;
#endif

typedef yakl::Array<real,yakl::memHost> realArrHost;

#if defined(__USE_HIP__)
  #include "hip/hip_runtime.h"
#endif

#include "params.h"

YAKL_INLINE real constexpr operator"" _fp( long double x ) {
  return static_cast<real>(x);
}

int constexpr ord      = 5;
int constexpr tord     = 3;
int constexpr hs       = (ord-1)/2;
int constexpr numState = 5;

int constexpr idR      = 0;  // "rho" (density)
int constexpr idU      = 1;  // "u" (u-wind)
int constexpr idV      = 2;  // "v" (v-wind)
int constexpr idW      = 3;  // "w" (w-wind)
int constexpr idT      = 4;  // thermodynamic variable
int constexpr idP      = 5;  // pressure

// Some physical constants
real constexpr PI    = 3.1415926535897932384626433832795028842;
real constexpr GRAV  = 9.8;
real constexpr CP    = 1004.;
real constexpr CV    = 717.;
real constexpr RD    = 287.;
real constexpr P0    = 1.0e5;
real constexpr C0    = 27.5629410929725921310572974482;
real constexpr GAMMA  = 1.40027894002789400278940027894;

template <class T> YAKL_INLINE T mymin( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> YAKL_INLINE T mymax( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

#endif
