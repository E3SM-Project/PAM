
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


YAKL_INLINE void diffTransformEulerX( SArray<real,numState,tord,tord> &state, 
                                      SArray<real,numState,tord,tord> &flux ,
                                      real pbar,
                                      SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord,tord> tmp_ruu, tmp_ruv, tmp_ruw, tmp_up, tmp_u;
  real tot_ruu, tot_ruv, tot_ruw, tot_up, tot_u;

  // Zero out intermediate arrays
  tmp_ruu = 0;
  tmp_ruv = 0;
  tmp_ruw = 0;
  tmp_up  = 0;
  tmp_u   = 0;

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real r = state(idR,0,ii);
    real u = state(idU,0,ii) / r;
    real v = state(idV,0,ii) / r;
    real w = state(idW,0,ii) / r;
    real p = state(idP,0,ii);

    tmp_ruu(0,ii) = r*u*u;
    tmp_ruv(0,ii) = r*u*v;
    tmp_ruw(0,ii) = r*u*w;
    tmp_up (0,ii) =   u*p;
    tmp_u  (0,ii) =   u;

    flux(idR,0,ii) = r*u;
    flux(idU,0,ii) = r*u*u + p;
    flux(idV,0,ii) = r*u*v;
    flux(idW,0,ii) = r*u*w;
    flux(idP,0,ii) =   u*p + (GAMMA-1)*pbar*u;
  }

  // Loop over the time derivatives
  for (int kt=0; kt<tord-1; kt++) {
    // Compute the state at the next time level
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        real d_dx = 0;
        for (int s=0; s<tord; s++) {
          d_dx += deriv(s,ii) * flux(l,kt,s);
        }
        state(l,kt+1,ii) = -d_dx/(kt+1._fp);
      }
    }

    // Compute ru* at the next time level
    for (int ii=0; ii<tord; ii++) {
      tot_ruu = 0;
      tot_ruv = 0;
      tot_ruw = 0;
      tot_up  = 0;
      tot_u   = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_ruu += state(idU,rt,ii) * state(idU,kt+1-rt,ii) - state(idR,rt,ii) * tmp_ruu(kt+1-rt,ii);
        tot_ruv += state(idU,rt,ii) * state(idV,kt+1-rt,ii) - state(idR,rt,ii) * tmp_ruv(kt+1-rt,ii);
        tot_ruw += state(idU,rt,ii) * state(idW,kt+1-rt,ii) - state(idR,rt,ii) * tmp_ruw(kt+1-rt,ii);
        tot_up  += state(idU,rt,ii) * state(idP,kt+1-rt,ii) - state(idR,rt,ii) * tmp_up (kt+1-rt,ii);
        tot_u   += state(idU,rt,ii)                         - state(idR,rt,ii) * tmp_u  (kt+1-rt,ii);
      }
      tmp_ruu(kt+1,ii) = tot_ruu / state(idR,0,ii);
      tmp_ruv(kt+1,ii) = tot_ruv / state(idR,0,ii);
      tmp_ruw(kt+1,ii) = tot_ruw / state(idR,0,ii);
      tmp_up (kt+1,ii) = tot_up  / state(idR,0,ii);
      tmp_u  (kt+1,ii) = tot_u   / state(idR,0,ii);

      // Compute the fluxes at the next time level
      flux(idR,kt+1,ii) = state(idU,kt+1,ii);
      flux(idU,kt+1,ii) = tmp_ruu(kt+1,ii) + state(idP,kt+1,ii);
      flux(idV,kt+1,ii) = tmp_ruv(kt+1,ii);
      flux(idW,kt+1,ii) = tmp_ruw(kt+1,ii);
      flux(idP,kt+1,ii) = tmp_up (kt+1,ii) + (GAMMA-1)*pbar*tmp_u(kt+1,ii);
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
