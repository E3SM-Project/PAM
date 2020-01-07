
#ifndef __ABCOEFS_H__
#define __ABCOEFS_H__

#include "grid.h"

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
void abcoefs(real1d const &dt3, int nstep, int na, int nb, int nc, real &at, real &bt, real &ct);

#endif

