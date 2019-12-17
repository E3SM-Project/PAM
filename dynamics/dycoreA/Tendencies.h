
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Riemann.h"
#include "Domain.h"
#include "Exchange.h"

class Tendencies {

  realArr stateLimits;
  realArr flux_r;
  realArr flux_re;
  realArr dph;
  SArray<real,tord> gllWts;
  SArray<real,ord,tord> to_gll;
  SArray<real,ord,tord> to_derivX_gll;
  SArray<real,ord,tord> to_derivY_gll;
  SArray<real,ord,tord> to_derivZ_gll;
  SArray<real,tord,tord> aderDerivX;
  SArray<real,tord,tord> aderDerivY;
  SArray<real,tord,tord> aderDerivZ;
  SArray<real,ord,ord,ord> wenoRecon;
  SArray<real,hs+2> wenoIdl;
  real wenoSigma;

public :


  void initialize(Domain const &dom);

  // Transform ord stencil cell averages into tord GLL point values
  YAKL_INLINE void reconStencil(SArray<real,ord> const &stencil, SArray<real,tord> &gll, int const doWeno,
                                SArray<real,ord,ord,ord> const &wenoRecon, SArray<real,ord,tord> const &to_gll,
                                SArray<real,hs+2> const &wenoIdl, real wenoSigma) {
    SArray<real,ord> coefs;
    if (doWeno) {
      compute_weno_coefs(wenoRecon,stencil,coefs,wenoIdl,wenoSigma);
    } else {
      for (int ii=0; ii<ord; ii++) {
        coefs(ii) = stencil(ii);
      }
    }

    for (int ii=0; ii<tord; ii++) {
      gll(ii) = 0.;
      for (int s=0; s<ord; s++) {
        gll(ii) += to_gll(s,ii) * coefs(s);
      }
    }
  }

  void compEulerTend_X(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

  void compEulerTend_Y(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

  void compEulerTend_Z(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

  void compEulerTend_S(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

  void stateBoundariesZ(realArr &state, Domain const &dom);

  void edgeBoundariesZ(realArr &stateLimits, realArr &fluxLimits, Domain const &dom);

  void compStrakaTend(realArr &state, Domain const &dom, Exchange &exch, Parallel const &par, realArr &tend);

};

#endif

