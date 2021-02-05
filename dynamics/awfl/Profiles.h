
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


  YAKL_INLINE real init_supercell_temperature(real z, real z_0, real z_trop, real z_top,
                                                      real T_0, real T_trop, real T_top) {
    if (z <= z_trop) {
      real lapse = - (T_trop - T_0) / (z_trop - z_0);
      return T_0 - lapse * (z - z_0);
    } else {
      real lapse = - (T_top - T_trop) / (z_top - z_trop);
      return T_trop - lapse * (z - z_trop);
    }
  }


  YAKL_INLINE real init_supercell_pressure_dry(real z, real z_0, real z_trop, real z_top,
                                                       real T_0, real T_trop, real T_top,
                                                       real p_0, real R_d) {
    if (z <= z_trop) {
      real lapse = - (T_trop - T_0) / (z_trop - z_0);
      real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
      return p_0 * pow( T / T_0 , GRAV/(R_d*lapse) );
    } else {
      // Get pressure at the tropopause
      real lapse = - (T_trop - T_0) / (z_trop - z_0);
      real p_trop = p_0 * pow( T_trop / T_0 , GRAV/(R_d*lapse) );
      // Get pressure at requested height
      lapse = - (T_top - T_trop) / (z_top - z_trop);
      if (lapse != 0) {
        real T = init_supercell_temperature(z, z_0, z_trop, z_top, T_0, T_trop, T_top);
        return p_trop * pow( T / T_trop , GRAV/(R_d*lapse) );
      } else {
        return p_trop * exp(-GRAV*(z-z_trop)/(R_d*T_trop));
      }
    }
  }

  
  YAKL_INLINE real init_supercell_relhum(real z, real z_0, real z_trop) {
    if (z <= z_trop) {
      return 1._fp - 0.75_fp * pow(z / z_trop , 1.25_fp );
    } else {
      return 0.25_fp;
    }
  }

  
  YAKL_INLINE real init_supercell_relhum_d_dz(real z, real z_0, real z_trop) {
    if (z <= z_trop) {
      return -0.9375_fp*pow(z/z_trop, 0.25_fp)/z_trop;
    } else {
      return 0;
    }
  }


  YAKL_INLINE real init_supercell_sat_mix_dry( real press , real T ) {
    return 380/(press) * exp( 17.27_fp * (T-273)/(T-36) );
  }


  YAKL_INLINE real init_supercell_sat_mix_dry_d_dT( real p , real T ) {
    return 2205033.6*exp(17.27*T/(T - 36) - 6424.44/(T - 36))/(p*(T*T - 72*T + 1296));
  }


  YAKL_INLINE real init_supercell_sat_mix_dry_d_dp( real p , real T ) {
    return -380*exp(17.27*T/(T - 36) - 6424.44/(T - 36))/(p*p);
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

