
#include "Kokkos_Core.hpp"
#include "abcoefs.h"

extern real at, bt, ct;
extern int  na, nb, nc, nstep;
extern real dt3[3];    // This cannot be defined as real *dt3

typedef Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > realArr;

// Compute the coefficients for the Adams-Bashforth scheme
extern "C" void abcoefs() {
  // Wrap pointers in unmanaged Kokkos Views
  real *dt3_r = dt3;
  realArr dt3(dt3_r,3);

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

