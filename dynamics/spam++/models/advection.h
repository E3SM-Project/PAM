#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "model.h"

// ********* Advection-Specific Stuff *************** //

// Compile-time constants
int constexpr nqdofs = 1;
int data_init_cond = DATA_INIT_GAUSSIAN;
int wind_init_cond = WIND_INIT_UNIFORM_X;

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
