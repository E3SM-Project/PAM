
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


YAKL_INLINE void diffTransformEulerConsX( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux, SArray<real,tord,tord> const &deriv ) {
  SArray<real,tord,tord> ruu, ruv, ruw, rut, rtgamma;
  real tot_ruu, tot_ruv, tot_ruw, tot_rut, tot_rtgamma;

  // Zero out intermediate arrays
  ruu     = 0;
  ruv     = 0;
  ruw     = 0;
  rut     = 0;
  rtgamma = 0;

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real r = state(idR ,0,ii);
    real u = state(idRU,0,ii) / r;
    real v = state(idRV,0,ii) / r;
    real w = state(idRW,0,ii) / r;
    real t = state(idRT,0,ii) / r;

    ruu    (0,ii) = r*u*u;
    ruv    (0,ii) = r*u*v;
    ruw    (0,ii) = r*u*w;
    rut    (0,ii) = r*u*t;
    rtgamma(0,ii) = pow( r*t , GAMMA );

    flux(idR ,0,ii) = r*u;
    flux(idRU,0,ii) = ruu(0,ii) + C0*rtgamma(0,ii);
    flux(idRV,0,ii) = ruv(0,ii);
    flux(idRW,0,ii) = ruw(0,ii);
    flux(idRT,0,ii) = rut(0,ii);
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
      tot_rut = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_ruu += state(idRU,rt,ii) * state(idRU,kt+1-rt,ii) - state(idR,rt,ii) * ruu(kt+1-rt,ii);
        tot_ruv += state(idRU,rt,ii) * state(idRV,kt+1-rt,ii) - state(idR,rt,ii) * ruv(kt+1-rt,ii);
        tot_ruw += state(idRU,rt,ii) * state(idRW,kt+1-rt,ii) - state(idR,rt,ii) * ruw(kt+1-rt,ii);
        tot_rut += state(idRU,rt,ii) * state(idRT,kt+1-rt,ii) - state(idR,rt,ii) * rut(kt+1-rt,ii);
      }
      ruu(kt+1,ii) = tot_ruu / state(idR,0,ii);
      ruv(kt+1,ii) = tot_ruv / state(idR,0,ii);
      ruw(kt+1,ii) = tot_ruw / state(idR,0,ii);
      rut(kt+1,ii) = tot_rut / state(idR,0,ii);

      // Compute rtgamma at the next time level
      tot_rtgamma = 0;
      for (int rt=0; rt<=kt; rt++) {
        tot_rtgamma += (kt+1._fp -rt) * ( GAMMA*rtgamma(rt,ii)*state(idRT,kt+1-rt,ii) - state(idRT,rt,ii)*rtgamma(kt+1-rt,ii) );
      }
      rtgamma(kt+1,ii) = ( GAMMA*rtgamma(0,ii)*state(idRT,kt+1,ii) + tot_rtgamma / (kt+1._fp) ) / state(idRT,0,ii);

      // Compute the fluxes at the next time level
      flux(idR ,kt+1,ii) = state(idRU,kt+1,ii);
      flux(idRU,kt+1,ii) = ruu(kt+1,ii) + C0*rtgamma(kt+1,ii)/2;
      flux(idRV,kt+1,ii) = ruv(kt+1,ii);
      flux(idRW,kt+1,ii) = ruw(kt+1,ii);
      flux(idRT,kt+1,ii) = rut(kt+1,ii);
    }
  }
}


YAKL_INLINE void diffTransformEulerConsY( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux, SArray<real,tord,tord> const &deriv ) {
  SArray<real,tord,tord> rvu, rvv, rvw, rvt, rtgamma;
  real tot_rvu, tot_rvv, tot_rvw, tot_rvt, tot_rtgamma;

  // Zero out intermediate arrays
  rvu     = 0;
  rvv     = 0;
  rvw     = 0;
  rvt     = 0;
  rtgamma = 0;

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real r = state(idR ,0,ii);
    real u = state(idRU,0,ii) / r;
    real v = state(idRV,0,ii) / r;
    real w = state(idRW,0,ii) / r;
    real t = state(idRT,0,ii) / r;

    rvu    (0,ii) = r*v*u;
    rvv    (0,ii) = r*v*v;
    rvw    (0,ii) = r*v*w;
    rvt    (0,ii) = r*v*t;
    rtgamma(0,ii) = pow( r*t , GAMMA );

    flux(idR ,0,ii) = r*v;
    flux(idRU,0,ii) = rvu(0,ii);
    flux(idRV,0,ii) = rvv(0,ii) + C0*rtgamma(0,ii);
    flux(idRW,0,ii) = rvw(0,ii);
    flux(idRT,0,ii) = rvt(0,ii);
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

    // Compute rv* at the next time level
    for (int ii=0; ii<tord; ii++) {
      tot_rvu = 0;
      tot_rvv = 0;
      tot_rvw = 0;
      tot_rvt = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_rvu += state(idRV,rt,ii) * state(idRU,kt+1-rt,ii) - state(idR,rt,ii) * rvu(kt+1-rt,ii);
        tot_rvv += state(idRV,rt,ii) * state(idRV,kt+1-rt,ii) - state(idR,rt,ii) * rvv(kt+1-rt,ii);
        tot_rvw += state(idRV,rt,ii) * state(idRW,kt+1-rt,ii) - state(idR,rt,ii) * rvw(kt+1-rt,ii);
        tot_rvt += state(idRV,rt,ii) * state(idRT,kt+1-rt,ii) - state(idR,rt,ii) * rvt(kt+1-rt,ii);
      }
      rvu(kt+1,ii) = tot_rvu / state(idR,0,ii);
      rvv(kt+1,ii) = tot_rvv / state(idR,0,ii);
      rvw(kt+1,ii) = tot_rvw / state(idR,0,ii);
      rvt(kt+1,ii) = tot_rvt / state(idR,0,ii);

      // Compute rtgamma at the next time level
      tot_rtgamma = 0;
      for (int rt=0; rt<=kt; rt++) {
        tot_rtgamma += (kt+1._fp -rt) * ( GAMMA*rtgamma(rt,ii)*state(idRT,kt+1-rt,ii) - state(idRT,rt,ii)*rtgamma(kt+1-rt,ii) );
      }
      rtgamma(kt+1,ii) = ( GAMMA*rtgamma(0,ii)*state(idRT,kt+1,ii) + tot_rtgamma / (kt+1._fp) ) / state(idRT,0,ii);

      // Compute the fluxes at the next time level
      flux(idR ,kt+1,ii) = state(idRV,kt+1,ii);
      flux(idRU,kt+1,ii) = rvu(kt+1,ii);
      flux(idRV,kt+1,ii) = rvv(kt+1,ii) + C0*rtgamma(kt+1,ii)/2;
      flux(idRW,kt+1,ii) = rvw(kt+1,ii);
      flux(idRT,kt+1,ii) = rvt(kt+1,ii);
    }
  }
}


YAKL_INLINE void diffTransformEulerConsZ( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux,
                                 SArray<real,tord,tord> &source, SArray<real,tord,tord> const &deriv, SArray<real,tord> const &hyRHOT, SArray<real,tord> const &hyRHO ) {
  SArray<real,tord,tord> rwu, rwv, rww, rwt, rtgamma;
  real tot_rwu, tot_rwv, tot_rww, tot_rwt, tot_rtgamma;

  // Zero out intermediate arrays
  rwu     = 0;
  rwv     = 0;
  rww     = 0;
  rwt     = 0;
  rtgamma = 0;

  // Compute the zeroth-order DTs of the intermediate functions and fluxes
  for (int ii=0; ii<tord; ii++) {
    real r = state(idR ,0,ii);
    real u = state(idRU,0,ii) / r;
    real v = state(idRV,0,ii) / r;
    real w = state(idRW,0,ii) / r;
    real t = state(idRT,0,ii) / r;

    rwu    (0,ii) = r*w*u;
    rwv    (0,ii) = r*w*v;
    rww    (0,ii) = r*w*w;
    rwt    (0,ii) = r*w*t;
    rtgamma(0,ii) = pow( r*t , GAMMA );

    flux(idR ,0,ii) = r*w;
    flux(idRU,0,ii) = rwu(0,ii);
    flux(idRV,0,ii) = rwv(0,ii);
    flux(idRW,0,ii) = rww(0,ii) + C0*rtgamma(0,ii) - C0*pow(hyRHOT(ii),GAMMA);
    flux(idRT,0,ii) = rwt(0,ii);

    source(0,ii) = -( state(idR,0,ii) - hyRHO(ii) )*GRAV;
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
        state(l,kt+1,ii) = -d_dx/(kt+1);
      }
    }
    for (int ii=0; ii<tord; ii++) {
      state(idRW,kt+1,ii) += source(kt,ii)/(kt+1);;
    }

    // Compute rw* at the next time level
    for (int ii=0; ii<tord; ii++) {
      tot_rwu = 0;
      tot_rwv = 0;
      tot_rww = 0;
      tot_rwt = 0;
      for (int rt=0; rt<=kt+1; rt++) {
        tot_rwu += state(idRW,rt,ii) * state(idRU,kt+1-rt,ii) - state(idR,rt,ii) * rwu(kt+1-rt,ii);
        tot_rwv += state(idRW,rt,ii) * state(idRV,kt+1-rt,ii) - state(idR,rt,ii) * rwv(kt+1-rt,ii);
        tot_rww += state(idRW,rt,ii) * state(idRW,kt+1-rt,ii) - state(idR,rt,ii) * rww(kt+1-rt,ii);
        tot_rwt += state(idRW,rt,ii) * state(idRT,kt+1-rt,ii) - state(idR,rt,ii) * rwt(kt+1-rt,ii);
      }
      rwu(kt+1,ii) = tot_rwu / state(idR,0,ii);
      rwv(kt+1,ii) = tot_rwv / state(idR,0,ii);
      rww(kt+1,ii) = tot_rww / state(idR,0,ii);
      rwt(kt+1,ii) = tot_rwt / state(idR,0,ii);

      // Compute rtgamma at the next time level
      tot_rtgamma = 0;
      for (int rt=0; rt<=kt; rt++) {
        tot_rtgamma += (kt+1-rt) * ( GAMMA*rtgamma(rt,ii)*state(idRT,kt+1-rt,ii) - state(idRT,rt,ii)*rtgamma(kt+1-rt,ii) );
      }
      rtgamma(kt+1,ii) = ( GAMMA*rtgamma(0,ii)*state(idRT,kt+1,ii) + tot_rtgamma / (kt+1) ) / state(idRT,0,ii);

      // Compute the fluxes at the next time level
      flux(idR ,kt+1,ii) = state(idRW,kt+1,ii);
      flux(idRU,kt+1,ii) = rwu(kt+1,ii);
      flux(idRV,kt+1,ii) = rwv(kt+1,ii);
      flux(idRW,kt+1,ii) = rww(kt+1,ii) + C0*rtgamma(kt+1,ii)/2;
      flux(idRT,kt+1,ii) = rwt(kt+1,ii);
      source(kt+1,ii) = -state(idR,kt+1,ii)*GRAV;
    }
  }
}


YAKL_INLINE void timeAvg( SArray<real,numState,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int l=0; l<numState; l++) {
      for (int ii=0; ii<tord; ii++) {
        dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1._fp);
      }
    }
    dtmult *= dom.dt;
  }
}


YAKL_INLINE void timeAvg( SArray<real,tord,tord> &dts , Domain const &dom ) {
  real dtmult = dom.dt;
  for (int kt=1; kt<tord; kt++) {
    for (int ii=0; ii<tord; ii++) {
      dts(0,ii) += dts(kt,ii) * dtmult / (kt+1._fp);
    }
    dtmult *= dom.dt;
  }
}


#endif
