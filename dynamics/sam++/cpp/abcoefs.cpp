
#include "abcoefs.h"

// Compute the Adams-Bashforth coefficients, at, bt, and ct
// 
// The na, nb, and nc indices are shifted to continue using the same
// three memory slots without having to copy data from one slot to another.
// It is typical practice for multi-step methods.
// na, nb, nc can be: {1,2,3}, {3,1,2}, or {2,3,1}
// dt3 contains the last three time step sizes
// 
// Boot-strapping procedure
// On the first time step, do forward Euler (no other data to use)
// On the second time step, do 2nd-order A-B
// On the third and continuing time steps, do 3rd-order A-B
void abcoefs(Domain &dom) {
  // coefficients for the Adams-Bashforth scheme
  real alpha, beta;

  if (dom.nstep >= 3) {          // Third time step and onward
    alpha = dom.dt3(dom.nb-1) / dom.dt3(dom.na-1);
    beta  = dom.dt3(dom.nc-1) / dom.dt3(dom.na-1);
    dom.ct = (2.+3.* alpha) / (6.* (alpha + beta) * beta);
    dom.bt = -(1.+2.*(alpha + beta) * dom.ct)/(2. * alpha);
    dom.at = 1. - dom.bt - dom.ct;
  } else if (dom.nstep >= 2) {   // Second time step
    dom.at = 3./2.;
    dom.bt = -1./2.;
    dom.ct = 0.;
  } else {                       // First time step
    dom.at = 1.;
    dom.bt = 0.;
    dom.ct = 0.;
  }
}


