#pragma once

#include <cmath>
#include <iostream>
//#include <cstring>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include "mpi.h"
#include <math.h>
#include "yaml-cpp/yaml.h"

////////////// These control the settings for SPAM++    //////////////

// Spatial derivatives order of accuracy ie Hodge stars [2,4,6] (vert only supports 2 for now)
uint constexpr diff_ord = 2;
uint constexpr vert_diff_ord = 2;

// Reconstruction types and order
enum class RECONSTRUCTION_TYPE { CFV, WENO, WENOFUNC };

RECONSTRUCTION_TYPE constexpr reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr dual_reconstruction_type = RECONSTRUCTION_TYPE::WENOFUNC;
uint constexpr dual_reconstruction_order = 5;

RECONSTRUCTION_TYPE constexpr coriolis_reconstruction_type = RECONSTRUCTION_TYPE::CFV;
uint constexpr coriolis_reconstruction_order = 3;

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
uint constexpr ic_quad_pts = 5;

// FIX THIS
// Halo sizes
uint constexpr maxhalosize = 15; //mymax(reconstruction_order+1,differential_order)/2; // IS THIS ALWAYS CORRECT?
uint constexpr mirroringhalo = 9; //mymax(reconstruction_order+1,differential_order)/2; // IS THIS ALWAYS CORRECT?

//0 = RKSimple, 1=SSPRK
#define _TIME_TYPE 0

//////////////////////////////////////////////////////////////////////////////////////////////

typedef unsigned long ulong;
typedef unsigned int  uint;

// Declaring the precision for the model
typedef double real;
#define REAL_NC NC_DOUBLE
#define REAL_MPI MPI_DOUBLE

#define PNETCDF_PUT_VAR ncmpi_put_vara_double
#define PNETCDF_PUT_VAR_ALL ncmpi_put_vara_double_all


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


#ifdef _LAYER
#include "layermodel-header.h"
#elif _EXTRUDED
#include "extrudedmodel-header.h"
#elif _ADVECTION
#include "advectionmodel-header.h"
#endif

template <uint nvars> void set_dofs_arr(SArray<int,2, nvars, 3> &dofs_arr, int var, int basedof, int extdof, int ndofs)
{
  dofs_arr(var, 0) = basedof;
  dofs_arr(var, 1) = extdof;
  dofs_arr(var, 2) = ndofs;
}



