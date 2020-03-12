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

// Declaring the precision for the model
typedef double real;
#define PNETCDF_PUT_VAR ncmpi_put_vara_double
#define PNETCDF_PUT_VAR_ALL ncmpi_put_vara_double_all
#define REAL_NC NC_DOUBLE
#define REAL_MPI MPI_DOUBLE


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



enum class RECONSTRUCTION_TYPE { CFV, UFV, WENO };
enum class TIME_TYPE { KGRK, ADER };
enum class GEOM_TYPE { UNIFORM_RECT, DISTORTED };

class Parameters
{
public:

  int nx_glob = -1;
  int ny_glob = -1;
  int nz_glob = -1;

  int Nsteps = -1;
  int Nout = -1;
  real dt = -1.;
  int Nstat = -1;
  std::string outputName = "output.nc";

  real etime;
//THESE ARE REALLY SPECIFIC TO UNIFORM RECT GEOM...
  real xlen, ylen, zlen;
  real xc, yc, zc;

};






// Number of Dimensions
uint constexpr number_of_dims = 2;

// Spatial order of accuracy for the model
uint constexpr differential_order = 8;

// Reconstruction type
RECONSTRUCTION_TYPE constexpr reconstruction_type = RECONSTRUCTION_TYPE::UFV;
uint constexpr reconstruction_order = 1;

// Halo sizes
uint maxhalosize = mymax(reconstruction_order,differential_order)/2; // IS THIS ALWAYS CORRECT?

// initial condition quadrature pts
uint constexpr ic_quad_pts = 3;

// Time scheme
TIME_TYPE constexpr time_type = TIME_TYPE::KGRK;
uint constexpr n_time_stages = 4;

// Grid geometry
GEOM_TYPE constexpr geom_type = GEOM_TYPE::UNIFORM_RECT;

#include "model-compile-consts.h"









#endif
