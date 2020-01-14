
#include "abcoefs.h"

// Compute the coefficients for the Adams-Bashforth scheme
extern "C" void abcoefs(int na, int nb, int nc, int nstep, real &dt3_p, real &at, real &bt, real &ct) {
  // Wrap pointers in unmanaged Kokkos Views
  umgReal1d dt3(&dt3_p,3);

  if (nstep >= 3) {
    real alpha = dt3(nb-1) / dt3(na-1);
    real beta  = dt3(nc-1) / dt3(na-1);
    ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta);
    bt = -(1.+2.*(alpha + beta) * ct)/(2. * alpha);
    at = 1. - bt - ct;
  } else if (nstep >= 2) {
    at = 3./2.;
    bt = -1./2.;
    ct = 0.;
  } else {
    at = 1.;
    bt = 0.;
    ct = 0.;
  }
}

