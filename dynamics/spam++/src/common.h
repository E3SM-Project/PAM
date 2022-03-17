#pragma once

#include <cmath>
#include <iostream>
#include <cstring>

#include "YAKL.h"
//#include <thread>

using yakl::SArray;
using yakl::c::parallel_for;
using yakl::c::SimpleBounds;
using yakl::c::Bounds;
using yakl::memDevice;
using yakl::memHost;
using yakl::styleC;
using yakl::index_t;

typedef unsigned long ulong;
typedef unsigned int  uint;

// Declaring the precision for the model
typedef double real;

#define PNETCDF_PUT_VAR ncmpi_put_vara_double
#define PNETCDF_PUT_VAR_ALL ncmpi_put_vara_double_all
#define REAL_NC NC_DOUBLE
#define REAL_MPI MPI_DOUBLE


// Define array types, presumed on the "device"
#if defined(__USE_CUDA__) || defined(__USE_HIP__)

typedef yakl::Array<real,1,yakl::memDevice,yakl::styleC> real1d;
typedef yakl::Array<real,2,yakl::memDevice,yakl::styleC> real2d;
typedef yakl::Array<real,3,yakl::memDevice,yakl::styleC> real3d;
typedef yakl::Array<real,4,yakl::memDevice,yakl::styleC> real4d;
typedef yakl::Array<real,5,yakl::memDevice,yakl::styleC> real5d;
typedef yakl::Array<real,6,yakl::memDevice,yakl::styleC> real6d;
typedef yakl::Array<real,7,yakl::memDevice,yakl::styleC> real7d;
typedef yakl::Array<real,8,yakl::memDevice,yakl::styleC> real8d;

typedef yakl::Array<int,1,yakl::memDevice,yakl::styleC> int1d;
typedef yakl::Array<int,2,yakl::memDevice,yakl::styleC> int2d;
typedef yakl::Array<int,3,yakl::memDevice,yakl::styleC> int3d;
typedef yakl::Array<int,4,yakl::memDevice,yakl::styleC> int4d;
typedef yakl::Array<int,5,yakl::memDevice,yakl::styleC> int5d;
typedef yakl::Array<int,6,yakl::memDevice,yakl::styleC> int6d;
typedef yakl::Array<int,7,yakl::memDevice,yakl::styleC> int7d;
typedef yakl::Array<int,8,yakl::memDevice,yakl::styleC> int8d;

typedef yakl::Array<bool,1,yakl::memDevice,yakl::styleC> bool1d;
typedef yakl::Array<bool,2,yakl::memDevice,yakl::styleC> bool2d;
typedef yakl::Array<bool,3,yakl::memDevice,yakl::styleC> bool3d;
typedef yakl::Array<bool,4,yakl::memDevice,yakl::styleC> bool4d;
typedef yakl::Array<bool,5,yakl::memDevice,yakl::styleC> bool5d;
typedef yakl::Array<bool,6,yakl::memDevice,yakl::styleC> bool6d;
typedef yakl::Array<bool,7,yakl::memDevice,yakl::styleC> bool7d;
typedef yakl::Array<bool,8,yakl::memDevice,yakl::styleC> bool8d;

#else

typedef yakl::Array<real,1,yakl::memHost,yakl::styleC> real1d;
typedef yakl::Array<real,2,yakl::memHost,yakl::styleC> real2d;
typedef yakl::Array<real,3,yakl::memHost,yakl::styleC> real3d;
typedef yakl::Array<real,4,yakl::memHost,yakl::styleC> real4d;
typedef yakl::Array<real,5,yakl::memHost,yakl::styleC> real5d;
typedef yakl::Array<real,6,yakl::memHost,yakl::styleC> real6d;
typedef yakl::Array<real,7,yakl::memHost,yakl::styleC> real7d;
typedef yakl::Array<real,8,yakl::memHost,yakl::styleC> real8d;

typedef yakl::Array<int,1,yakl::memHost,yakl::styleC> int1d;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleC> int2d;
typedef yakl::Array<int,3,yakl::memHost,yakl::styleC> int3d;
typedef yakl::Array<int,4,yakl::memHost,yakl::styleC> int4d;
typedef yakl::Array<int,5,yakl::memHost,yakl::styleC> int5d;
typedef yakl::Array<int,6,yakl::memHost,yakl::styleC> int6d;
typedef yakl::Array<int,7,yakl::memHost,yakl::styleC> int7d;
typedef yakl::Array<int,8,yakl::memHost,yakl::styleC> int8d;

typedef yakl::Array<bool,1,yakl::memHost,yakl::styleC> bool1d;
typedef yakl::Array<bool,2,yakl::memHost,yakl::styleC> bool2d;
typedef yakl::Array<bool,3,yakl::memHost,yakl::styleC> bool3d;
typedef yakl::Array<bool,4,yakl::memHost,yakl::styleC> bool4d;
typedef yakl::Array<bool,5,yakl::memHost,yakl::styleC> bool5d;
typedef yakl::Array<bool,6,yakl::memHost,yakl::styleC> bool6d;
typedef yakl::Array<bool,7,yakl::memHost,yakl::styleC> bool7d;
typedef yakl::Array<bool,8,yakl::memHost,yakl::styleC> bool8d;

#endif

// Define array types on the "host"


typedef yakl::Array<real,1,yakl::memHost,yakl::styleC> realHost1d;
typedef yakl::Array<real,2,yakl::memHost,yakl::styleC> realHost2d;
typedef yakl::Array<real,3,yakl::memHost,yakl::styleC> realHost3d;
typedef yakl::Array<real,4,yakl::memHost,yakl::styleC> realHost4d;
typedef yakl::Array<real,5,yakl::memHost,yakl::styleC> realHost5d;
typedef yakl::Array<real,6,yakl::memHost,yakl::styleC> realHost6d;
typedef yakl::Array<real,7,yakl::memHost,yakl::styleC> realHost7d;
typedef yakl::Array<real,8,yakl::memHost,yakl::styleC> realHost8d;

typedef yakl::Array<int,1,yakl::memHost,yakl::styleC> intHost1d;
typedef yakl::Array<int,2,yakl::memHost,yakl::styleC> intHost2d;
typedef yakl::Array<int,3,yakl::memHost,yakl::styleC> intHost3d;
typedef yakl::Array<int,4,yakl::memHost,yakl::styleC> intHost4d;
typedef yakl::Array<int,5,yakl::memHost,yakl::styleC> intHost5d;
typedef yakl::Array<int,6,yakl::memHost,yakl::styleC> intHost6d;
typedef yakl::Array<int,7,yakl::memHost,yakl::styleC> intHost7d;
typedef yakl::Array<int,8,yakl::memHost,yakl::styleC> intHost8d;

typedef yakl::Array<bool,1,yakl::memHost,yakl::styleC> boolHost1d;
typedef yakl::Array<bool,2,yakl::memHost,yakl::styleC> boolHost2d;
typedef yakl::Array<bool,3,yakl::memHost,yakl::styleC> boolHost3d;
typedef yakl::Array<bool,4,yakl::memHost,yakl::styleC> boolHost4d;
typedef yakl::Array<bool,5,yakl::memHost,yakl::styleC> boolHost5d;
typedef yakl::Array<bool,6,yakl::memHost,yakl::styleC> boolHost6d;
typedef yakl::Array<bool,7,yakl::memHost,yakl::styleC> boolHost7d;
typedef yakl::Array<bool,8,yakl::memHost,yakl::styleC> boolHost8d;

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
uint constexpr vert_diff_ord = _VERT_DIFF_ORDER;

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

uint constexpr mirroringhalo = 9; //mymax(reconstruction_order+1,differential_order)/2; // IS THIS ALWAYS CORRECT?

#include "model-compile-consts.h"
