
#pragma once

#include "phys_params.h"

namespace profiles {

  YAKL_INLINE real initConstTheta_density(real t0, real z, real Rd, real cp, real gamma, real p0, real C0) {
    real exner = 1._fp - GRAV*z/(cp*t0);
    real p = pow( exner , cp/Rd ) * p0;
    real rt = pow( p/C0 , 1._fp/gamma );
    return rt / t0;
  }


  YAKL_INLINE real initConstTheta_pressure(real t0, real z, real Rd, real cp, real gamma, real p0, real C0) {
    real r = initConstTheta_density(t0,z,Rd,cp,gamma,p0,C0);
    return C0*pow(r*t0,gamma);
  }


  YAKL_INLINE real initConstTheta_pressureDeriv(real t0, real z, real Rd, real cp, real gamma, real p0, real C0) {
    real p = initConstTheta_pressure(t0,z,Rd,cp,gamma,p0,C0);
    return -GRAV/(t0*Rd)*pow(p0,Rd/cp)*pow(p,-Rd/cp+1);
  }


  /*
    Gives a linear ellipsiod centered at (x0,y0,z0) with radius (xrad,yrad,zrad) and amplitude amp
  */
  YAKL_INLINE real ellipsoid_linear(real x   , real y   , real z ,
                                    real x0  , real y0  , real z0,
                                    real xrad, real yrad, real zrad, real amp) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = sqrt( xn*xn + yn*yn + zn*zn );
    return amp * max( 1._fp - dist , 0._fp );
  }


  /*
    Gives a cosine ellipsiod centered at (x0,y0,z0) with radius (xrad,yrad,zrad) and amplitude amp
  */
  YAKL_INLINE real ellipsoid_cosine(real x   , real y   , real z ,
                                    real x0  , real y0  , real z0,
                                    real xrad, real yrad, real zrad, real amp, real pwr) {
    real val = 0;
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real zn = (z-z0)/zrad;
    real dist = sqrt( xn*xn + yn*yn + zn*zn );
    if (dist <= 1._fp) {
      val = amp * pow( (cos(M_PI*dist)+1)/2 , pwr );
    }
    return val;
  }

}

