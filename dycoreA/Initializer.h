
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Hydrostasis.h"
#include "Exchange.h"
#include "TransformMatrices.h"
#include "TimeIntegrator.h"
#include "mpi.h"
#include "Domain.h"
#include "Parallel.h"


void initialize_mpi( int *argc , char ***argv , Parallel &par );

void initialize(realArr &state, Domain &dom, Parallel &par, Exchange &exch, TimeIntegrator &tint);

YAKL_INLINE real ellipsoid_linear(real const x   , real const y   , real const z ,
                                  real const x0  , real const y0  , real const z0,
                                  real const xrad, real const yrad, real const zrad, real const amp) {
  real xn = (x-x0)/xrad;
  real yn = (y-y0)/yrad;
  real zn = (z-z0)/zrad;
  real dist = sqrt( xn*xn + yn*yn + zn*zn );
  return amp * mymax( 1._fp - dist , 0._fp );
}


YAKL_INLINE real ellipsoid_cosine(real const x   , real const y   , real const z ,
                                  real const x0  , real const y0  , real const z0,
                                  real const xrad, real const yrad, real const zrad, real const amp, real const pwr) {
  real val = 0;
  real xn = (x-x0)/xrad;
  real yn = (y-y0)/yrad;
  real zn = (z-z0)/zrad;
  real dist = sqrt( xn*xn + yn*yn + zn*zn );
  if (dist <= 1._fp) {
    val = amp * pow( (cos(PI*dist)+1)/2 , pwr );
  }
  return val;
}


#endif
