
#ifndef _CONST_H_
#define _CONST_H_

#include <cmath>
#include "Array.h"
#include "SArray.h"
#if defined(__USE_HIP__)
  #include "hip/hip_runtime.h"
#endif

using yakl::SArray;

// Declaring the precision for the model
typedef float         real;

// The parameters depend on the real type
#include "params.h"

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

// Spatial order of accuracy for the model
int constexpr ord      = 5;

// Time order of accuracy for the model
int constexpr tord     = 3;

// The number of halo cells needed for stencils
int constexpr hs       = (ord-1)/2;

// The number of entries in the state vector
int constexpr numState = 5;

// Locations of different variables within the state vector. Pressure is added in the
// vertical direction for reconstructed state values and derivatives, for the time DTs, 
// and for the stateLimits array (and the subsequent flux difference splitting). This is
// Required in the vertical direction to maintain hydrostasis in the flux difference
// splitting. 
int constexpr idR      = 0;  // "rho" (density)
int constexpr idU      = 1;  // "u" (u-wind)
int constexpr idV      = 2;  // "v" (v-wind)
int constexpr idW      = 3;  // "w" (w-wind)
int constexpr idT      = 4;  // thermodynamic variable
int constexpr idP      = 5;  // pressure

// Some physical constants
real constexpr PI    = 3.1415926535897932384626433832795028842; // obvious
real constexpr GRAV  = 9.8;                                     // acceleration due to gravity
real constexpr CP    = 1004.;                                   // specific heat at constant pressure
real constexpr CV    = 717.;                                    // specific heat at constant volume
real constexpr RD    = 287.;                                    // dry air constant
real constexpr P0    = 1.0e5;                                   // pressure for normalized Exner pressure
real constexpr C0    = 27.5629410929725921310572974482;         // constant of eqn. of state p=C0*pow(rho*theta,GAMMA)
real constexpr GAMMA  = 1.40027894002789400278940027894;        // CP/CV

// Specifying templated min and max functions
template <class T> YAKL_INLINE T mymin( T const v1 , T const v2 ) {
  if (v1 < v2) { return v1; }
  else         { return v2; }
}
template <class T> YAKL_INLINE T mymax( T const v1 , T const v2 ) {
  if (v1 > v2) { return v1; }
  else         { return v2; }
}

#endif
