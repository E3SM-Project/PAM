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


// Add mode for various operators
enum class ADD_MODE { REPLACE, ADD };



class Parameters
{
public:

  int nx_glob = -1;
  int ny_glob = -1;
  int nz_glob = -1;

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

};




enum class RECONSTRUCTION_TYPE { CFV, WENO, WENOFUNC };
enum class TIME_TYPE { KGRK, SSPRK };
enum class GEOM_TYPE { UNIFORM_RECT, };


// SHOULD SET ALL OF THIS AT COMPILE TIME BASED ON THE MODEL VIA COMMAND LINE ARGUMENTS/FLAGS WITH DEFAULTS!
// OR IT SHOULD REALLY MOVE INTO THE model-compile-const.h ?
// ideally run scripts acquire basically all of these, via -D...!

// what sorts of models do we have initially
// advection 1D/2D/3D (does both primal and dual advection, including implied PV advection with various choices for Q...)
// layeredmodel1D
// layeredmodel2D- CAN THIS BE MERGED WITH 1D?
// splitmodel2D
// splitmodel3D- CAN THIS BE MERGED WITH 2D?

// eventually add support for out of slice velocity/etc. in layered models and splitmodel2D


// Number of Dimensions
uint constexpr ndims = 2;

//THIS SHOULD BE HODGE STAR ORDER!
//MAYBE EVEN SET DIFFERENT STARS WITH DIFFERENT ACCURACIES?
// Spatial order of accuracy for the model
uint constexpr diff_ord = 2;

//EVENTUALLY DISTINGUISH BETWEEN DENSITY AND DENSITYFCT RECONS HERE...

// Reconstruction type
RECONSTRUCTION_TYPE constexpr reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr coriolis_reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr coriolis_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr vert_reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr vert_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_vert_reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_vert_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr coriolis_vert_reconstruction_type = RECONSTRUCTION_TYPE::CFV;
uint constexpr coriolis_vert_reconstruction_order = 3;

// How to handle PV flux term
// ADD AL81-TYPE SCHEME HERE EVENTUALLY AS WELL
enum class QF_MODE { EC, NOEC };
QF_MODE constexpr qf_choice = QF_MODE::EC;

// initial condition quadrature pts
uint constexpr ic_quad_pts = 3;

// Time scheme
TIME_TYPE constexpr time_type = TIME_TYPE::SSPRK;
uint constexpr n_time_stages = 3;

// Grid geometry
GEOM_TYPE constexpr geom_type = GEOM_TYPE::UNIFORM_RECT;


// FIX THIS
// Halo sizes
uint constexpr maxhalosize = 15; //mymax(reconstruction_order+1,differential_order)/2; // IS THIS ALWAYS CORRECT?



#include "model-compile-consts.h"




#endif
