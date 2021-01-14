
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


  YAKL_INLINE real init_supercell_theta(real zloc, real trop_a, real trop_b, real trop_c,
                                        real strat_a, real strat_b, real z_tr) {
    if (zloc <= z_tr) {
      return exp(trop_a*zloc*zloc + trop_b*zloc + trop_c);
    } else {
      return exp(strat_a*zloc + strat_b);
    }
  }


  YAKL_INLINE real init_supercell_exner_trop(real zloc , real trop_a, real trop_b, real trop_c, real cp) {
    return 1 - 1.0/2.0*sqrt(M_PI)*(-GRAV*exp((1.0/4.0)*pow(trop_b, 2)/trop_a)*erf((1.0/2.0)*trop_b/sqrt(trop_a)) + GRAV*exp((1.0/4.0)*pow(trop_b, 2)/trop_a)*erf((1.0/2.0)*(2*trop_a*zloc + trop_b)/sqrt(trop_a)))*exp(-trop_c)/(cp*sqrt(trop_a));
  }


  YAKL_INLINE real init_supercell_exner_strat(real zloc, real strat_a, real strat_b, real cp, real z_tr) {
    return (GRAV*exp(strat_a*z_tr) - GRAV*exp(strat_a*zloc))*exp(-strat_a*z_tr - strat_a*zloc - strat_b)/(cp*strat_a);
  }


  YAKL_INLINE real init_supercell_exner(real zloc, real trop_a, real trop_b, real trop_c,
                                        real strat_a, real strat_b, real cp, real z_tr) {
    if (zloc <= z_tr) {
      return init_supercell_exner_trop(zloc, trop_a, trop_b, trop_c, cp);
    } else {
      return init_supercell_exner_trop(z_tr, trop_a, trop_b, trop_c, cp) + 
             init_supercell_exner_strat(zloc, strat_a, strat_b, cp, z_tr);
    }
  }


  // YAKL_INLINE real init_supercell_pressure(real zloc, real trop_a, real trop_b, real trop_c,
  //                                          real strat_a, real strat_b, real cp, real Rd, real p0) {
  //   real p = pow( exner , cp/Rd ) * p0;
  //   return p;
  // }


  // YAKL_INLINE real init_supercell_density(real zloc, real a, real b, real c, real cp, real Rd, real p0, real C0, real gamma) {
  //   real p = init_supercell_pressure( zloc , a , b , c , cp , Rd , p0 );
  //   real theta = init_supercell_theta( zloc , a , b , c );
  //   return pow( p / C0 , 1/gamma ) / theta;
  // }


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

