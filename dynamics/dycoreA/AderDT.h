
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


// Computes time derivatives of the state, state derivatves, mass flux, and energy flux
// Dimensions of state, derive, rflux, and reflux are (time DT order, spatial GLL point)
// state    : tord time DTs of the state at tord GLL points across a cell
// deriv    : tord time DTs of the state spatial derivative at tord GLL points across a cell
// rflux    : tord time DTs of the mass flux (rho*u) at tord GLL points across a cell
// reflux   : tord time DTs of the energy flux (rho*e*u+u*p) at tord GLL points across a cell
// utend    : tord time DTs of u*du/dx in the DTs and then -u*du/dx-(dp/dx)/rho at the function exit
// vtend    : tord time DTs of u*dv/dx in the DTs and then -u*dv/dx             at the function exit
// wtend    : tord time DTs of u*dw/dx in the DTs and then -u*dw/dx             at the function exit
// deriv_mat: Matrix that transforms tord GLL points into tord GLL points of the spatial derivative
YAKL_INLINE void diffTransformEulerX( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real         ,tord,tord> &utend,
                                      SArray<real         ,tord,tord> &vtend,
                                      SArray<real         ,tord,tord> &wtend,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> rflux;     // r*u
  SArray<real,tord,tord> tmp_u_dp;  // u*dp/dx
  SArray<real,tord,tord> tmp_rr_dp; // (1/rho)*dp/dx
  SArray<real,tord,tord> tmp_p_du;  // p*du/dx
  SArray<real,tord,tord> tmp_re;    // rho*e
  SArray<real,tord,tord> tmp_u_re;  // u*rho*e
  SArray<real,tord,tord> tmp_u_p;   // u*p

  // The term (dp/dx)/rho involves division. Because of this, the time DTs
  // must be zeroed out to kill the term that would otherwise look like recursion
  for (int kt=1; kt<tord; kt++) {
    for int (ii=0; ii<tord; ii++) {
      tmp_rr_dp(kt,ii) = 0;
    }
  }

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    // state values
    real r  = state(idR,0,ii);
    real u  = state(idU,0,ii);
    real v  = state(idV,0,ii);
    real w  = state(idW,0,ii);
    real p  = state(idT,0,ii);
    // state derivatives
    real dr = deriv(idR,0,ii);
    real du = deriv(idU,0,ii);
    real dv = deriv(idV,0,ii);
    real dw = deriv(idW,0,ii);
    real dp = deriv(idT,0,ii);

    real ke = 0.5_fp*r*(u*u+v*v+w+w);  // kinetic energy
    real re = p*CV/RD+ke;              // rho*e
    
    // Initialize the 0th-order time DTs (i.e., the values)
    rflux    (0,ii) = r*u;
    utend    (0,ii) = u*du;
    vtend    (0,ii) = u*dv;
    wtend    (0,ii) = u*dw;
    tmp_u_dp (0,ii) = u*dp;
    tmp_rr_dp(0,ii) = dp/r;
    tmp_p_du (0,ii) = p*du;
    tmp_re   (0,ii) = re;
    tmp_u_re (0,ii) = u*re;
    tmp_u_p  (0,ii) = u*p;
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and p
    state(idU,kt+1,ii) = -(utend   (kt,ii) +      tmp_rr_dp(kt,ii))/(kt+1);
    state(idV,kt+1,ii) = -(vtend   (kt,ii)                        )/(kt+1);
    state(idW,kt+1,ii) = -(wtend   (kt,ii)                        )/(kt+1);
    state(idP,kt+1,ii) = -(tmp_u_dp(kt,ii) + GAMMA*tmp_p_du(kt,ii))/(kt+1);

    // Compute (kt+1)th DT of rho and rho*e
    for (int ii=0; ii<tord; ii++) {
      real drflux_dx  = 0;
      real dreflux_dx = 0;
      // Matrix-vector multiply against the spatial differentiation matrix
      for (int s=0; s<tord; s++) {
        drflux_dx  += deriv(s,ii) * rflux (kt,s);
        dreflux_dx += deriv(s,ii) * ( tmp_u_re(kt,ii) + tmp_u_p(kt,ii) );
      }
      state(idR,kt+1,ii) = -drflux_dx /(kt+1);
      tmp_re(   kt+1,ii) = -dreflux_dx/(kt+1);
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dx = 0;
        // Matrix-vector multiply against the spatial differentiation matrix
        for (int s=0; s<tord; s++) {
          d_dx += deriv(s,ii) * state(l,kt,s);
        }
        deriv(l,kt+1,ii) = d_dx;
      }
    }

    // Compute the (kt+1)th DT of all temporary variables
    // Nearly all of these are of the form f*g
    // Except for (dp/dx)/rho, which is of the form f/g
    for (int ii=0; ii<tord; ii++) {
      real tot_rflux     = 0;
      real tot_tmp_u_du  = 0;
      real tot_tmp_u_dv  = 0;
      real tot_tmp_u_dw  = 0;
      real tot_tmp_u_dp  = 0;
      real tot_tmp_p_du  = 0;
      real tot_tmp_u_re  = 0;
      real tot_tmp_u_p   = 0;
      real tot_tmp_rr_dp = 0;
      for (int rt=0; tr<=kt+1; rt++) {
        tot_rflux     += state(idR,rt,ii) * state(idU,kt+1-rt,ii);
        tot_tmp_u_du  += state(idU,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_u_dv  += state(idU,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_u_dw  += state(idU,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_u_dp  += state(idU,rt,ii) * deriv(idP,kt+1-rt,ii);
        tot_tmp_p_du  += state(idP,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_u_re  += state(idU,rt,ii) * tmp_re(   kt+1-rt,ii);
        tot_tmp_u_p   += state(idU,rt,ii) * state(idP,kt+1-rt,ii);
        tot_tmp_rr_dp += deriv(idP,rt,ii) - state(idR,rt,ii) * tmp_rr_dp(kt+1-rt,ii);
      }
      rflux    (kt+1,ii) = tot_rflux   ;
      utend    (kt+1,ii) = tot_tmp_u_du;
      vtend    (kt+1,ii) = tot_tmp_u_dv;
      wtend    (kt+1,ii) = tot_tmp_u_dw;
      tmp_u_dp (kt+1,ii) = tot_tmp_u_dp;
      tmp_p_du (kt+1,ii) = tot_tmp_p_du;
      tmp_u_re (kt+1,ii) = tot_tmp_u_re;
      tmp_u_p  (kt+1,ii) = tot_tmp_u_p ;
      tmp_rr_dp(kt+1,ii) = tot_tmp_rr_dp / state(idR,0,ii);
    }
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      utend(kt,ii) = -utend(kt,ii) - rr_dp;
      vtend(kt,ii) = -vtend(kt,ii);
      wtend(kt,ii) = -wtend(kt,ii);
    }
  }
}



YAKL_INLINE void diffTransformEulerY( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &flux ,
                                      real pbar,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_rvu, tmp_rvv, tmp_rvw, tmp_vp, tmp_v;
  real tot_rvu, tot_rvv, tot_rvw, tot_vp, tot_v;

  // Zero out intermediate arrays
  tmp_rvu = 0;
  tmp_rvv = 0;
  tmp_rvw = 0;
  tmp_vp  = 0;
  tmp_v   = 0;

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real r = state(idR,0,ii);
    real u = state(idU,0,ii) / r;
    real v = state(idV,0,ii) / r;
    real w = state(idW,0,ii) / r;
    real p = state(idP,0,ii);

    tmp_rvu(0,ii) = r*v*u;
    tmp_rvv(0,ii) = r*v*v;
    tmp_rvw(0,ii) = r*v*w;
    tmp_vp (0,ii) =   v*p;
    tmp_v  (0,ii) =   v;

    flux(idR,0,ii) = r*v;
    flux(idU,0,ii) = r*v*u;
    flux(idV,0,ii) = r*v*v + p;
    flux(idW,0,ii) = r*v*w;
    flux(idP,0,ii) =   v*p + (GAMMA-1)*pbar*v;
  }

  // Loop over the time derivatives
  for (int kt=0; kt<tord-1; kt++) {
    // Compute the state at the next time level
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dy = 0;
        for (int s=0; s<tord; s++) {
          d_dy += deriv(s,ii) * flux(l,kt,s);
        }
        state(l,kt+1,ii) = -d_dy/(kt+1._fp);
      }
    }

    // Compute rv* at the next time level
    for (int ii=0; ii<tord; ii++) {
      tot_rvu = 0;
      tot_rvv = 0;
      tot_rvw = 0;
      tot_vp  = 0;
      tot_v   = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_rvu += state(idV,rt,ii) * state(idU,kt+1-rt,ii) - state(idR,rt,ii) * tmp_rvu(kt+1-rt,ii);
        tot_rvv += state(idV,rt,ii) * state(idV,kt+1-rt,ii) - state(idR,rt,ii) * tmp_rvv(kt+1-rt,ii);
        tot_rvw += state(idV,rt,ii) * state(idW,kt+1-rt,ii) - state(idR,rt,ii) * tmp_rvw(kt+1-rt,ii);
        tot_vp  += state(idV,rt,ii) * state(idP,kt+1-rt,ii) - state(idR,rt,ii) * tmp_vp (kt+1-rt,ii);
        tot_v   += state(idV,rt,ii)                         - state(idR,rt,ii) * tmp_v  (kt+1-rt,ii);
      }
      tmp_rvu(kt+1,ii) = tot_rvu / state(idR,0,ii);
      tmp_rvv(kt+1,ii) = tot_rvv / state(idR,0,ii);
      tmp_rvw(kt+1,ii) = tot_rvw / state(idR,0,ii);
      tmp_vp (kt+1,ii) = tot_vp  / state(idR,0,ii);
      tmp_v  (kt+1,ii) = tot_v   / state(idR,0,ii);

      // Compute the fluxes at the next time level
      flux(idR,kt+1,ii) = state(idV,kt+1,ii);
      flux(idU,kt+1,ii) = tmp_rvu(kt+1,ii);
      flux(idV,kt+1,ii) = tmp_rvv(kt+1,ii) + state(idP,kt+1,ii);
      flux(idW,kt+1,ii) = tmp_rvw(kt+1,ii);
      flux(idP,kt+1,ii) = tmp_vp (kt+1,ii) + (GAMMA-1)*pbar*tmp_v(kt+1,ii);
    }
  }
}



YAKL_INLINE void diffTransformEulerZ( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &flux ,
                                      real pbar,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_rwu, tmp_rwv, tmp_rww, tmp_wp, tmp_w;
  real tot_rwu, tot_rwv, tot_rww, tot_wp, tot_w;

  // Zero out intermediate arrays
  tmp_rwu = 0;
  tmp_rwv = 0;
  tmp_rww = 0;
  tmp_wp  = 0;
  tmp_w   = 0;

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real r = state(idR,0,ii);
    real u = state(idU,0,ii) / r;
    real v = state(idV,0,ii) / r;
    real w = state(idW,0,ii) / r;
    real p = state(idP,0,ii);

    tmp_rwu(0,ii) = r*w*u;
    tmp_rwv(0,ii) = r*w*v;
    tmp_rww(0,ii) = r*w*w;
    tmp_wp (0,ii) =   w*p;
    tmp_w  (0,ii) =   w;

    flux(idR,0,ii) = r*w;
    flux(idU,0,ii) = r*w*u;
    flux(idV,0,ii) = r*w*v;
    flux(idW,0,ii) = r*w*w + p;
    flux(idP,0,ii) =   w*p + (GAMMA-1)*pbar*w;
  }

  // Loop over the time derivatives
  for (int kt=0; kt<tord-1; kt++) {
    // Compute the state at the next time level
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dy = 0;
        for (int s=0; s<tord; s++) {
          d_dy += deriv(s,ii) * flux(l,kt,s);
        }
        state(l,kt+1,ii) = -d_dy/(kt+1._fp);
      }
    }

    // Compute rv* at the next time level
    for (int ii=0; ii<tord; ii++) {
      tot_rwu = 0;
      tot_rwv = 0;
      tot_rww = 0;
      tot_wp  = 0;
      tot_w   = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_rwu += state(idW,rt,ii) * state(idU,kt+1-rt,ii) - state(idR,rt,ii) * tmp_rwu(kt+1-rt,ii);
        tot_rwv += state(idW,rt,ii) * state(idV,kt+1-rt,ii) - state(idR,rt,ii) * tmp_rwv(kt+1-rt,ii);
        tot_rww += state(idW,rt,ii) * state(idW,kt+1-rt,ii) - state(idR,rt,ii) * tmp_rww(kt+1-rt,ii);
        tot_wp  += state(idW,rt,ii) * state(idP,kt+1-rt,ii) - state(idR,rt,ii) * tmp_wp (kt+1-rt,ii);
        tot_w   += state(idW,rt,ii)                         - state(idR,rt,ii) * tmp_w  (kt+1-rt,ii);
      }
      tmp_rwu(kt+1,ii) = tot_rwu / state(idR,0,ii);
      tmp_rwv(kt+1,ii) = tot_rwv / state(idR,0,ii);
      tmp_rww(kt+1,ii) = tot_rww / state(idR,0,ii);
      tmp_wp (kt+1,ii) = tot_wp  / state(idR,0,ii);
      tmp_w  (kt+1,ii) = tot_w   / state(idR,0,ii);

      // Compute the fluxes at the next time level
      flux(idR,kt+1,ii) = state(idW,kt+1,ii);
      flux(idU,kt+1,ii) = tmp_rwu(kt+1,ii);
      flux(idV,kt+1,ii) = tmp_rwv(kt+1,ii);
      flux(idW,kt+1,ii) = tmp_rww(kt+1,ii) + state(idP,kt+1,ii);
      flux(idP,kt+1,ii) = tmp_wp (kt+1,ii) + (GAMMA-1)*pbar*tmp_w(kt+1,ii);
    }
  }
}




template <int N> YAKL_INLINE void timeAvg( SArray<real,N,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int l=0; l<N; l++) {
      for (int ii=0; ii<tord; ii++) {
        dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1);
      }
    }
    dtmult *= dom.dt;
  }
}


YAKL_INLINE void timeAvg( SArray<real,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      dts(0,ii) += dts(kt,ii) * dtmult / (kt+1);
    }
    dtmult *= dom.dt;
  }
}


#endif
