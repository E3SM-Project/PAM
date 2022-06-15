#pragma once

#include <cmath>
#include <iostream>
//#include <cstring>
#include "mpi.h"
#include "yaml-cpp/yaml.h"
#include <array>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>

using yakl::c::parallel_for;
using yakl::c::SimpleBounds;

typedef unsigned long ulong;
typedef unsigned int uint;

////////////// These control the settings for SPAM++    //////////////

// Declaring the precision for the model
typedef double real;
#define REAL_MPI MPI_DOUBLE
//#define REAL_NC NC_DOUBLE

// Spatial derivatives order of accuracy ie Hodge stars [2,4,6] (vert only
// supports 2 for now)
uint constexpr diff_ord = 2;
uint constexpr vert_diff_ord = 2;

// Reconstruction types and order
enum class RECONSTRUCTION_TYPE { CFV, WENO, WENOFUNC };

RECONSTRUCTION_TYPE constexpr reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr coriolis_reconstruction_type =
    RECONSTRUCTION_TYPE::CFV;
uint constexpr coriolis_reconstruction_order = 3;

RECONSTRUCTION_TYPE constexpr vert_reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr vert_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_vert_reconstruction_type =
    RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_vert_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr coriolis_vert_reconstruction_type =
    RECONSTRUCTION_TYPE::CFV;
uint constexpr coriolis_vert_reconstruction_order = 3;

// How to handle PV flux term
// ADD AL81-TYPE SCHEME HERE EVENTUALLY AS WELL
enum class QF_MODE { EC, NOEC };
QF_MODE constexpr qf_choice = QF_MODE::EC;

// initial condition quadrature pts
uint constexpr ic_quad_pts_x = 5;
uint constexpr ic_quad_pts_y = 5;
uint constexpr ic_quad_pts_z = 5;

// FIX THIS
// Halo sizes
uint constexpr maxhalosize =
    15; // mymax(reconstruction_order+1,differential_order)/2;
        // // IS THIS ALWAYS CORRECT?
uint constexpr mirroringhalo =
    9; // mymax(reconstruction_order+1,differential_order)/2;
       // // IS THIS ALWAYS CORRECT?

// 0 = RKSimple, 1=SSPRK
#define _TIME_TYPE 1

//////////////////////////////////////////////////////////////////////////////////////////////

real const pi =
    3.141592653589793238462643383279502884197169399375105820974944_fp;

//#define PNETCDF_PUT_VAR ncmpi_put_vara_double
//#define PNETCDF_PUT_VAR_ALL ncmpi_put_vara_double_all

// Specifying templated min and max functions
template <class T> YAKL_INLINE T mymin(T const v1, T const v2) {
  if (v1 < v2) {
    return v1;
  } else {
    return v2;
  }
}
template <class T> YAKL_INLINE T mymax(T const v1, T const v2) {
  if (v1 > v2) {
    return v1;
  } else {
    return v2;
  }
}

// Utility function for settings dofs_arr
template <uint nvars>
void set_dofs_arr(SArray<int, 2, nvars, 3> &dofs_arr, int var, int basedof,
                  int extdof, int ndofs) {
  dofs_arr(var, 0) = basedof;
  dofs_arr(var, 1) = extdof;
  dofs_arr(var, 2) = ndofs;
}

// Utility functions for serial output
void inline debug_print(std::string out, int masterproc) {
#ifdef PAM_DEBUG
  if (masterproc) {
    std::cout << out << "\n";
  }
#endif
}
void inline serial_print(std::string out, int masterproc) {
  if (masterproc) {
    std::cout << out << "\n";
  }
}

class Parameters {
public:
  int nx_glob = -1;
  int ny_glob = -1;
  int nz_dual = -1;
  int nens = -1;

  int Nsteps = -1;
  int Nout = -1;
  real dtcrm = -1.;
  real dtphys = -1.;
  int crm_per_phys = -1;
  int Nstat = -1;
  std::string outputName;
  std::string tstype;

  real xlen, ylen, zlen;
  real xc, yc, zc;

  int masterproc;
};

// Add mode for various operators
enum class ADD_MODE { REPLACE, ADD };

// Boundary types
enum class BND_TYPE { PERIODIC, NONE };

#if defined _HAMILTONIAN && defined _LAYER
#include "layermodel-common.h"
#elif defined _HAMILTONIAN && defined _EXTRUDED
#include "extrudedmodel-common.h"
#elif defined _ADVECTION && defined _LAYER
#include "layeradvection-common.h"
#elif defined _ADVECTION && defined _EXTRUDED
#include "extrudedadvection-common.h"
#endif
