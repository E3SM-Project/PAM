
#ifndef _RIEMANN_H_
#define _RIEMANN_H_

#include "const.h"
#include "SArray.h"


YAKL_INLINE void godunovLinearX(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                SArray<real,numState> &f1, SArray<real,numState> &f2,
                                SArray<real,numState> &upw ) {
  SArray<real,numState> ch1, ch2, chu, ev;
  real const windTol = 0; // 1.0e-6;

  // Compute interface values
  real r = 0.5_fp * ( s1(idR) + s2(idR) );
  real u = 0.5_fp * ( s1(idU) + s2(idU) ) / r;
  real v = 0.5_fp * ( s1(idV) + s2(idV) ) / r;
  real w = 0.5_fp * ( s1(idW) + s2(idW) ) / r;
  real p = 0.5_fp * ( s1(idP) + s2(idP) );
  real cs = sqrt( GAMMA * p / r );

  // Compute left and right characteristic variables
  ch1(0) =  u*f1(0)/(2*cs) - f1(1)/(2*cs) + f1(4)/(2*t);
  ch1(1) = -u*f1(0)/(2*cs) + f1(1)/(2*cs) + f1(4)/(2*t);
  ch1(2) = f1(0) -   f1(4)/t;
  ch1(3) = f1(2) - v*f1(4)/t;
  ch1(4) = f1(3) - w*f1(4)/t;

  ch2(0) =  u*f2(0)/(2*cs) - f2(1)/(2*cs) + f2(4)/(2*t);
  ch2(1) = -u*f2(0)/(2*cs) + f2(1)/(2*cs) + f2(4)/(2*t);
  ch2(2) = f2(0) -   f2(4)/t;
  ch2(3) = f2(2) - v*f2(4)/t;
  ch2(4) = f2(3) - w*f2(4)/t;

  ev(0) = u-cs;
  ev(1) = u+cs;
  ev(2) = u;
  ev(3) = u;
  ev(4) = u;

  // Compute the upwind characteristics
  for (int l=0; l<numState; l++) {
    if        (ev(l) >= windTol) {
      chu(l) = ch1(l);
    } else if (ev(l) < windTol) {
      chu(l) = ch2(l);
    } else {
      chu(l) = 0.5_fp * (ch1(l) + ch2(l));
    }
  }

  // Compute the fluxes
  upw(0) = chu(0) + chu(1) + chu(2);
  upw(1) = (u-cs)*chu(0) + (u+cs)*chu(1) + u*chu(2);
  upw(2) = v*chu(0) + v*chu(1) + chu(3);
  upw(3) = w*chu(0) + w*chu(1) + chu(4);
  upw(4) = t*chu(0) + t*chu(1);
}


YAKL_INLINE void godunovLinearY(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                SArray<real,numState> &f1, SArray<real,numState> &f2,
                                SArray<real,numState> &upw ) {
  SArray<real,numState> ch1, ch2, chu, ev;
  real const windTol = 0; // 1.0e-6;

  // Compute interface values
  real r = 0.5_fp * ( s1(idR ) + s2(idR ) );
  real u = 0.5_fp * ( s1(idRU) + s2(idRU) ) / r;
  real v = 0.5_fp * ( s1(idRV) + s2(idRV) ) / r;
  real w = 0.5_fp * ( s1(idRW) + s2(idRW) ) / r;
  real t = 0.5_fp * ( s1(idRT) + s2(idRT) ) / r;
  real p = C0 * pow( r*t , GAMMA );
  real cs = sqrt( GAMMA * p / r );

  // Compute left and right characteristic variables
  ch1(0) =  v*f1(0)/(2*cs) - f1(2)/(2*cs) + f1(4)/(2*t);
  ch1(1) = -v*f1(0)/(2*cs) + f1(2)/(2*cs) + f1(4)/(2*t);
  ch1(2) = f1(0) -   f1(4)/t;
  ch1(3) = f1(1) - u*f1(4)/t;
  ch1(4) = f1(3) - w*f1(4)/t;

  ch2(0) =  v*f2(0)/(2*cs) - f2(2)/(2*cs) + f2(4)/(2*t);
  ch2(1) = -v*f2(0)/(2*cs) + f2(2)/(2*cs) + f2(4)/(2*t);
  ch2(2) = f2(0) -   f2(4)/t;
  ch2(3) = f2(1) - u*f2(4)/t;
  ch2(4) = f2(3) - w*f2(4)/t;

  ev(0) = v-cs;
  ev(1) = v+cs;
  ev(2) = v;
  ev(3) = v;
  ev(4) = v;

  // Compute the upwind characteristics
  for (int l=0; l<numState; l++) {
    if        (ev(l) >= windTol) {
      chu(l) = ch1(l);
    } else if (ev(l) < windTol) {
      chu(l) = ch2(l);
    } else {
      chu(l) = 0.5_fp * (ch1(l) + ch2(l));
    }
  }

  // Compute the fluxes
  upw(0) = chu(0) + chu(1) + chu(2);
  upw(1) = u*chu(0) + u*chu(1) + chu(3);
  upw(2) = (v-cs)*chu(0) + (v+cs)*chu(1) + v*chu(2);
  upw(3) = w*chu(0) + w*chu(1) + chu(4);
  upw(4) = t*chu(0) + t*chu(1);
}


YAKL_INLINE void godunovLinearZ(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                SArray<real,numState> &f1, SArray<real,numState> &f2,
                                SArray<real,numState> &upw ) {
  SArray<real,numState> ch1, ch2, chu, ev;
  real const windTol = 0; // 1.0e-6;

  // Compute interface values
  real r = 0.5_fp * ( s1(idR ) + s2(idR ) );
  real u = 0.5_fp * ( s1(idRU) + s2(idRU) ) / r;
  real v = 0.5_fp * ( s1(idRV) + s2(idRV) ) / r;
  real w = 0.5_fp * ( s1(idRW) + s2(idRW) ) / r;
  real t = 0.5_fp * ( s1(idRT) + s2(idRT) ) / r;
  real p = C0 * pow( r*t , GAMMA );
  real cs = sqrt( GAMMA * p / r );

  // Compute left and right characteristic variables
  ch1(0) = f1(0) - f1(4)/t;
  ch1(1) = f1(1) - f1(4)*u/t;
  ch1(2) = f1(2) - f1(4)*v/t;
  ch1(3) =  f1(0)*w/(2*cs) - f1(3)/(2*cs) + f1(4)/(2*t);
  ch1(4) = -f1(0)*w/(2*cs) + f1(3)/(2*cs) + f1(4)/(2*t);

  ch2(0) = f2(0) - f2(4)/t;
  ch2(1) = f2(1) - f2(4)*u/t;
  ch2(2) = f2(2) - f2(4)*v/t;
  ch2(3) =  f2(0)*w/(2*cs) - f2(3)/(2*cs) + f2(4)/(2*t);
  ch2(4) = -f2(0)*w/(2*cs) + f2(3)/(2*cs) + f2(4)/(2*t);

  ev(0) = w;
  ev(1) = w;
  ev(2) = w;
  ev(3) = w-cs;
  ev(4) = w+cs;

  // Compute the upwind characteristics
  for (int l=0; l<numState; l++) {
    if        (ev(l) >= windTol) {
      chu(l) = ch1(l);
    } else if (ev(l) < windTol) {
      chu(l) = ch2(l);
    } else {
      chu(l) = 0.5_fp * (ch1(l) + ch2(l));
    }
  }

  // Compute the fluxes
  upw(0) = chu(0)   + chu(3)        + chu(4);
  upw(1) = chu(1)   + chu(3)*u      + chu(4)*u;
  upw(2) = chu(2)   + chu(3)*v      + chu(4)*v;
  upw(3) = chu(0)*w + chu(3)*(w-cs) + chu(4)*(w+cs);
  upw(4) =            chu(3)*t      + chu(4)*t;
}


YAKL_INLINE void riemannX(SArray<real,numState> &s1, SArray<real,numState> &s2,
                          SArray<real,numState> &f1, SArray<real,numState> &f2,
                          SArray<real,numState> &upw) {
  godunovLinearX(s1, s2, f1, f2, upw);
}


YAKL_INLINE void riemannY(SArray<real,numState> &s1, SArray<real,numState> &s2,
                          SArray<real,numState> &f1, SArray<real,numState> &f2,
                          SArray<real,numState> &upw) {
  godunovLinearY(s1, s2, f1, f2, upw);
}


YAKL_INLINE void riemannZ(SArray<real,numState> &s1, SArray<real,numState> &s2,
                          SArray<real,numState> &f1, SArray<real,numState> &f2,
                          SArray<real,numState> &upw) {
  godunovLinearZ(s1, s2, f1, f2, upw);
}


#endif
