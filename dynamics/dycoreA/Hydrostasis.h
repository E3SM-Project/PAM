
#ifndef _HYDROSTASIS_H_
#define _HYDROSTASIS_H_

#include "const.h"

namespace hydro {

  /*
    Computes the hydrostatic density (r) at the given constant
    potential temperature (t0) and location (zloc)
  */
  YAKL_INLINE void hydroConstTheta(real const t0, real const z, real &r) {
    real exner = 1._fp - GRAV*z/(CP*t0);
    real p = pow( exner , CP/RD ) * P0;
    real rt = pow( p/C0 , 1._fp/GAMMA );
    r = rt / t0;
  }

}

#endif
