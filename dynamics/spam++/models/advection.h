#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "driver.h"

// *********** COMPILE TIME CONSTANTS ************** //
// EVENTUALLY THESE SHOULD ALL BE COMPILE TIME/PRE-PROCESSOR FLAGS...
// WITH REASONABLE DEFAULTS!

// Number of Dimensions
int constexpr ndims = 1;

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
int constexpr maxhalosize_z = 0 // mymax(reconstruction_order_y,differential_order_y)/2; //IS THIS ALWAYS CORRECT?
int constexpr maxhalosize_y = 0 // mymax(reconstruction_order_z,differential_order_z)/2; //IS THIS ALWAYS CORRECT?

// initial condition quadrature pts
int constexpr ic_quad_pts_x = 3;
int constexpr ic_quad_pts_y = 1; //3;
int constexpr ic_quad_pts_z = 1; //3;

// Time scheme
int constexpr time_type = TIME_TYPE_KGRK;
int constexpr time_order = 2;

// Number of variables
int constexpr nprognostic = 1;
int constexpr nconstant = 1;
int constexpr ndiagnostic = 1;
int constexpr nqdofs = 1;

// Initial conditions
int data_init_cond = DATA_INIT_GAUSSIAN;
int wind_init_cond = WIND_INIT_UNIFORM_X;


// ************************* //

// Initial conditions related variables and functions

int constexpr DATA_INIT_GAUSSIAN   = 1;
int constexpr DATA_INIT_VORTICES   = 2;
int constexpr DATA_INIT_SQUARE     = 3;

int constexpr WIND_INIT_UNIFORM_X     = 1;
int constexpr WIND_INIT_UNIFORM_Y     = 2;
int constexpr WIND_INIT_UNIFORM_Z     = 3;
int constexpr WIND_INIT_DEFORMATIONAL = 4;

YAKL_INLINE real gaussian(real x, real y, real z) {
  return x;
}

YAKL_INLINE real vortices(real x, real y, real z) {
  return x;
}

YAKL_INLINE real square(real x, real y, real z) {
  return x;
}

YAKL_INLINE real uniform_x_wind(real x, real y, real z) {
  return x;
}

YAKL_INLINE real uniform_y_wind(real x, real y, real z) {
  return x;
}

YAKL_INLINE real uniform_z_wind(real x, real y, real z) {
  return x;
}

YAKL_INLINE real deformational_wind(real x, real y, real z) {
  return x;
}

#endif
