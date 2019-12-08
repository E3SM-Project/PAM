
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


YAKL_INLINE void diffTransformEulerX( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real         ,tord,tord> &utend,
                                      SArray<real         ,tord,tord> &vtend,
                                      SArray<real         ,tord,tord> &wtend,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_r_u;     // r*u
  SArray<real,tord,tord> tmp_u_dp;    // u*dp/dx
  SArray<real,tord,tord> tmp_rr_dp;   // (1/rho)*dp/dx
  SArray<real,tord,tord> tmp_p_du;    // p*du/dx

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
    
    // Initialize the 0th-order time DTs (i.e., the values)
    tmp_r_u  (0,ii) = r*u;
    utend    (0,ii) = u*du;
    vtend    (0,ii) = u*dv;
    wtend    (0,ii) = u*dw;
    tmp_u_dp (0,ii) = u*dp;
    tmp_rr_dp(0,ii) = dp/r;
    tmp_p_du (0,ii) = p*du;
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and p
    state(idU,kt+1,ii) = -(utend   (kt,ii) +      tmp_rr_dp(kt,ii))/(kt+1);  // u
    state(idV,kt+1,ii) = -(vtend   (kt,ii)                        )/(kt+1);  // v
    state(idW,kt+1,ii) = -(wtend   (kt,ii)                        )/(kt+1);  // w
    state(idT,kt+1,ii) = -(tmp_u_dp(kt,ii) + GAMMA*tmp_p_du(kt,ii))/(kt+1);  // p

    // Compute (kt+1)th DT of rho
    for (int ii=0; ii<tord; ii++) {
      real drflux_dx  = 0;
      // Matrix-vector multiply against the spatial differentiation matrix
      for (int s=0; s<tord; s++) {
        drflux_dx += deriv_mat(s,ii) * tmp_r_u(kt,s);
      }
      state(idR,kt+1,ii) = -drflux_dx/(kt+1);
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dx = 0;
        // Matrix-vector multiply against the spatial differentiation matrix
        for (int s=0; s<tord; s++) {
          d_dx += deriv_mat(s,ii) * state(l,kt,s);
        }
        deriv(l,kt+1,ii) = d_dx;
      }
    }

    // Compute the (kt+1)th DT of all temporary variables
    // Nearly all of these are of the form f*g
    // Except for (dp/dx)/rho, which is of the form f/g
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_r_u   = 0;
      real tot_tmp_u_du  = 0;
      real tot_tmp_u_dv  = 0;
      real tot_tmp_u_dw  = 0;
      real tot_tmp_u_dp  = 0;
      real tot_tmp_p_du  = 0;
      real tot_tmp_rr_dp = 0;
      for (int rt=0; tr<=kt+1; rt++) {
        tot_tmp_r_u   += state(idR,rt,ii) * state(idU,kt+1-rt,ii);
        tot_tmp_u_du  += state(idU,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_u_dv  += state(idU,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_u_dw  += state(idU,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_u_dp  += state(idU,rt,ii) * deriv(idP,kt+1-rt,ii);
        tot_tmp_p_du  += state(idP,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_rr_dp += deriv(idP,rt,ii) - state(idR,rt,ii) * tmp_rr_dp(kt+1-rt,ii);
      }
      tmp_r_u  (kt+1,ii) = tot_tmp_r_u ;
      utend    (kt+1,ii) = tot_tmp_u_du;
      vtend    (kt+1,ii) = tot_tmp_u_dv;
      wtend    (kt+1,ii) = tot_tmp_u_dw;
      tmp_u_dp (kt+1,ii) = tot_tmp_u_dp;
      tmp_p_du (kt+1,ii) = tot_tmp_p_du;
      tmp_rr_dp(kt+1,ii) = tot_tmp_rr_dp / state(idR,0,ii);
    } // ii-loop
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
                                      SArray<real,numState,tord,tord> &deriv,
                                      SArray<real         ,tord,tord> &utend,
                                      SArray<real         ,tord,tord> &vtend,
                                      SArray<real         ,tord,tord> &wtend,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_r_v;     // r*v
  SArray<real,tord,tord> tmp_v_dp;    // v*dp/dy
  SArray<real,tord,tord> tmp_rr_dp;   // (1/rho)*dp/dy
  SArray<real,tord,tord> tmp_p_dv;    // p*dv/dy

  // The term (dp/dy)/rho involves division. Because of this, the time DTs
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
    real du = deriv(idU,0,ii);
    real dv = deriv(idV,0,ii);
    real dw = deriv(idW,0,ii);
    real dp = deriv(idT,0,ii);
    
    // Initialize the 0th-order time DTs (i.e., the values)
    tmp_r_v  (0,ii) = r*v;
    utend    (0,ii) = v*du;
    vtend    (0,ii) = v*dv;
    wtend    (0,ii) = v*dw;
    tmp_v_dp (0,ii) = v*dp;
    tmp_rr_dp(0,ii) = dp/r;
    tmp_p_dv (0,ii) = p*dv;
  } // ii-loop

  // Loop over the time derivatives, computing the (kt+1)th time DTs in each iteration
  for (int kt=0; kt<tord-1; kt++) {

    // Compute (kt+1)th DT of u, v, w, and p
    state(idU,kt+1,ii) = -(utend   (kt,ii)                        )/(kt+1);  // u
    state(idV,kt+1,ii) = -(vtend   (kt,ii) +      tmp_rr_dp(kt,ii))/(kt+1);  // v
    state(idW,kt+1,ii) = -(wtend   (kt,ii)                        )/(kt+1);  // w
    state(idT,kt+1,ii) = -(tmp_v_dp(kt,ii) + GAMMA*tmp_p_dv(kt,ii))/(kt+1);  // p

    // Compute (kt+1)th DT of rho
    for (int ii=0; ii<tord; ii++) {
      real drflux_dy  = 0;
      // Matrix-vector multiply against the spatial differentiation matrix
      for (int s=0; s<tord; s++) {
        drflux_dy += deriv_mat(s,ii) * tmp_r_v(kt,s);
      }
      state(idR,kt+1,ii) = -drflux_dy/(kt+1);
    }

    // Compute spatial derivative of the (kt+1)th DTs of the state
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dy = 0;
        // Matrix-vector multiply against the spatial differentiation matrix
        for (int s=0; s<tord; s++) {
          d_dy += deriv_mat(s,ii) * state(l,kt,s);
        }
        deriv(l,kt+1,ii) = d_dy;
      }
    }

    // Compute the (kt+1)th DT of all temporary variables
    // Nearly all of these are of the form f*g
    // Except for (dp/dy)/rho, which is of the form f/g
    for (int ii=0; ii<tord; ii++) {
      real tot_tmp_r_v   = 0;
      real tot_tmp_v_du  = 0;
      real tot_tmp_v_dv  = 0;
      real tot_tmp_v_dw  = 0;
      real tot_tmp_v_dp  = 0;
      real tot_tmp_p_dv  = 0;
      real tot_tmp_rr_dp = 0;
      for (int rt=0; tr<=kt+1; rt++) {
        tot_tmp_r_v   += state(idR,rt,ii) * state(idV,kt+1-rt,ii);
        tot_tmp_v_du  += state(idV,rt,ii) * deriv(idU,kt+1-rt,ii);
        tot_tmp_v_dv  += state(idV,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_v_dw  += state(idV,rt,ii) * deriv(idW,kt+1-rt,ii);
        tot_tmp_v_dp  += state(idV,rt,ii) * deriv(idP,kt+1-rt,ii);
        tot_tmp_p_dv  += state(idP,rt,ii) * deriv(idV,kt+1-rt,ii);
        tot_tmp_rr_dp += deriv(idP,rt,ii) - state(idR,rt,ii) * tmp_rr_dp(kt+1-rt,ii);
      }
      tmp_r_v  (kt+1,ii) = tot_tmp_r_v ;
      utend    (kt+1,ii) = tot_tmp_v_du;
      vtend    (kt+1,ii) = tot_tmp_v_dv;
      wtend    (kt+1,ii) = tot_tmp_v_dw;
      tmp_v_dp (kt+1,ii) = tot_tmp_v_dp;
      tmp_p_dv (kt+1,ii) = tot_tmp_p_du;
      tmp_rr_dp(kt+1,ii) = tot_tmp_rr_dp / state(idR,0,ii);
    } // ii-loop
  } // kt-loop

  // Transform utend, vtend, and wtend into the RHS of the wind equations
  for (int kt=0; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      utend(kt,ii) = -utend(kt,ii);
      vtend(kt,ii) = -vtend(kt,ii) - rr_dp;
      wtend(kt,ii) = -wtend(kt,ii);
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
