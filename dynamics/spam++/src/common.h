#ifndef _COMMON_H_
#define _COMMON_H_

#include <cmath>
#include <iostream>
#include <cstring>

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


// Add mode for various operators
enum class ADD_MODE { REPLACE, ADD };

//Boundary types
enum class BND_TYPE { PERIODIC, NONE };



class Parameters
{
public:

  int nx_glob = -1;
  int ny_glob = -1;
  int nz = -1;

  int Nsteps = -1;
  int Nout = -1;
  real dt = -1.;
  real cfl = -1.;
  int Nstat = -1;
  std::string outputName = "output.nc";
  std::string TStype = "dt";

  real etime;
//THESE ARE REALLY SPECIFIC TO UNIFORM RECT GEOM...
  real xlen, ylen, zlen;
  real xc, yc, zc;

  std::string xbnd = "periodic";
  std::string ybnd = "periodic";

};






// These are all set in CMakeLists.txt by reading compile_consts.build

// Spatial order of accuracy for the model
uint constexpr diff_ord = _DIFF_ORDER;

// Reconstruction types and order
enum class RECONSTRUCTION_TYPE { CFV, WENO, WENOFUNC };

RECONSTRUCTION_TYPE constexpr reconstruction_type = RECONSTRUCTION_TYPE::_PRIMAL_RECON_TYPE;
uint constexpr reconstruction_order = _PRIMAL_RECON_ORDER;

RECONSTRUCTION_TYPE constexpr dual_reconstruction_type = RECONSTRUCTION_TYPE::_DUAL_RECON_TYPE;
uint constexpr dual_reconstruction_order = _DUAL_RECON_ORDER;

RECONSTRUCTION_TYPE constexpr coriolis_reconstruction_type = RECONSTRUCTION_TYPE::_CORIOLIS_RECON_TYPE;
uint constexpr coriolis_reconstruction_order = _CORIOLIS_RECON_ORDER;

RECONSTRUCTION_TYPE constexpr vert_reconstruction_type = RECONSTRUCTION_TYPE::_PRIMAL_VERT_RECON_TYPE;
uint constexpr vert_reconstruction_order = _PRIMAL_VERT_RECON_ORDER;

RECONSTRUCTION_TYPE constexpr dual_vert_reconstruction_type = RECONSTRUCTION_TYPE::_DUAL_VERT_RECON_TYPE;
uint constexpr dual_vert_reconstruction_order = _DUAL_VERT_RECON_ORDER;

RECONSTRUCTION_TYPE constexpr coriolis_vert_reconstruction_type = RECONSTRUCTION_TYPE::_CORIOLIS_VERT_RECON_TYPE;
uint constexpr coriolis_vert_reconstruction_order = _CORIOLIS_VERT_RECON_ORDER;

// How to handle PV flux term
// ADD AL81-TYPE SCHEME HERE EVENTUALLY AS WELL
enum class QF_MODE { EC, NOEC };
QF_MODE constexpr qf_choice = QF_MODE::_QF_CHOICE;

// initial condition quadrature pts
uint constexpr ic_quad_pts = _IC_QUAD_PTS;

// FIX THIS
// Halo sizes
uint constexpr maxhalosize = 15; //mymax(reconstruction_order+1,differential_order)/2; // IS THIS ALWAYS CORRECT?

#include "model-compile-consts.h"




#endif
