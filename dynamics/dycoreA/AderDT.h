
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


YAKL_INLINE void diffTransformEulerPrimX( SArray<real,numState,tord,tord> &state, 
                                          SArray<real,numState,tord,tord> &deriv,
                                          SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord> r, u, cs2ort, rt;  // rho, u, cs^2/(rho*theta), rho*theta
  // Precompute terms that will remain temporally constant ("Frozen Jacobian")
  for (int ii=0; ii<tord; ii++) {
    r     (ii) = state(idR ,0,ii);
    u     (ii) = state(idU ,0,ii);
    rt    (ii) = state(idRT,0,ii);
    real p   = C0*pow(rt(ii),GAMMA);
    real cs2 = GAMMA*p/r;
    cs2ort(ii) = cs2/(r*t);
  }
  // Loop over the time derivatives
  for (int kt=0; kt<tord-1; kt++) {
    // Compute the state at the next time level
    for (int ii=0; ii<tord; ii++) {
      // state time deriv      Slow advective terms      Fast non-advective terms
      state(idR ,kt+1,ii) = -( u(ii)*deriv(idR ,kt,ii) + r     (ii)*deriv(idU ,kt,ii) ) / (kt+1);
      state(idU ,kt+1,ii) = -( u(ii)*deriv(idU ,kt,ii) + cs2ort(ii)*deriv(idRT,kt,ii) ) / (kt+1);
      state(idV ,kt+1,ii) = -( u(ii)*deriv(idV ,kt,ii)                                ) / (kt+1);
      state(idW ,kt+1,ii) = -( u(ii)*deriv(idW ,kt,ii)                                ) / (kt+1);
      state(idRT,kt+1,ii) = -( u(ii)*deriv(idRT,kt,ii) + rt    (ii)*deriv(idU ,kt,ii) ) / (kt+1);
    }
    if (kt < tord-2) {
      // Comput the spatial derivative of the state at the next time level
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dx = 0;
          for (int s=0; s<tord; s++) {
            d_dx += deriv_mat(s,ii) * state(l,kt+1,s);
          }
          deriv(l,kt+1,ii) = -d_dx;
        }
      }
    }
  }
}



YAKL_INLINE void diffTransformEulerPrimY( SArray<real,numState,tord,tord> &state, 
                                          SArray<real,numState,tord> &deriv,
                                          SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord> r, v, cs2ort, rt;  // rho, u, cs^2/(rho*theta), rho*theta
  // Precompute terms that will remain temporally constant ("Frozen Jacobian")
  for (int ii=0; ii<tord; ii++) {
    r     (ii) = state(idR ,0,ii);
    v     (ii) = state(idV ,0,ii);
    rt    (ii) = state(idRT,0,ii);
    real p   = C0*pow(rt(ii),GAMMA);
    real cs2 = GAMMA*p/r;
    cs2ort(ii) = cs2/(r*t);
  }
  // Loop over the time derivatives
  for (int kt=0; kt<tord-1; kt++) {
    // Compute the state at the next time level
    for (int ii=0; ii<tord; ii++) {
      // state time deriv      Slow advective terms      Fast non-advective terms
      state(idR ,kt+1,ii) = -( v(ii)*deriv(idR ,kt,ii) + r     (ii)*deriv(idV ,kt,ii) ) / (kt+1);
      state(idU ,kt+1,ii) = -( v(ii)*deriv(idU ,kt,ii)                                ) / (kt+1);
      state(idV ,kt+1,ii) = -( v(ii)*deriv(idV ,kt,ii) + cs2ort(ii)*deriv(idRT,kt,ii) ) / (kt+1);
      state(idW ,kt+1,ii) = -( v(ii)*deriv(idW ,kt,ii)                                ) / (kt+1);
      state(idRT,kt+1,ii) = -( v(ii)*deriv(idRT,kt,ii) + rt    (ii)*deriv(idV ,kt,ii) ) / (kt+1);
    }
    if (kt < tord-2) {
      // Overwrite the derivative of the state
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dy = 0;
          for (int s=0; s<tord; s++) {
            d_dy += deriv_mat(s,ii) * state(l,kt+1,s);
          }
          deriv(l,ii) = -d_dy/(kt+1);
        }
      }
    }
  }
}



YAKL_INLINE void diffTransformEulerPrimZ( SArray<real,numState+1,tord,tord> &state, 
                                          SArray<real,numState+1,tord> &deriv,
                                          SArray<real,tord,tord> const &deriv_mat ) {
  SArray<real,tord> r, rr, w, cs2ot, rt;  // rho, 1/rho, u, cs^2/theta, rho*theta
  // Precompute terms that will remain temporally constant ("Frozen Jacobian")
  for (int ii=0; ii<tord; ii++) {
    r    (ii) = state(idR ,0,ii);
    w    (ii) = state(idW ,0,ii);
    rt   (ii) = state(idRT,0,ii);
    rr   (ii) = 1._fp/r(ii);
    real p   = C0*pow(rt(ii),GAMMA);
    real cs2 = GAMMA*p/r;
    cs2ot(ii) = w*cs2/t;
  }
  // Loop over the time derivatives
  for (int kt=0; kt<tord-1; kt++) {
    // Compute the state at the next time level
    for (int ii=0; ii<tord; ii++) {
      // state time deriv      Slow advective terms      Fast non-advective terms
      state(idR ,kt+1,ii) = -( w(ii)*deriv(idR ,kt,ii) + r (ii)*deriv(idW,kt,ii) ) / (kt+1);
      state(idU ,kt+1,ii) = -( w(ii)*deriv(idU ,kt,ii)                           ) / (kt+1);
      state(idV ,kt+1,ii) = -( w(ii)*deriv(idV ,kt,ii)                           ) / (kt+1);
      state(idW ,kt+1,ii) = -( w(ii)*deriv(idW ,kt,ii) + rr(ii)*deriv(idP,kt,ii) ) / (kt+1);
      state(idRT,kt+1,ii) = -( w(ii)*deriv(idRT,kt,ii) + rt(ii)*deriv(idW,kt,ii) ) / (kt+1);
      // Update pressure term as a scaling of the rho*theta update through the eqn of state
      state(idP ,kt+1,ii) = cs2ot(ii) * state(idRT,kt+1,ii);
    }
    if (kt < tord-2) {
      // Overwrite the derivative of the state
      for (int l=0; l<numState+1; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dz = 0;
          for (int s=0; s<tord; s++) {
            d_dz += deriv_mat(s,ii) * state(l,kt+1,s);
          }
          deriv(l,ii) = -d_dz/(kt+1);
        }
      }
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
