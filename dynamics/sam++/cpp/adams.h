
#ifndef _ADAMS_H_
#define _ADAMS_H_

#include "grid.h"

void adams(Domain const &dom, real1d const &dz, real2d const &rho, real2d const &rhow,
           real5d &dudt, real5d &dvdt, real5d dwdt, real4d &u, real4d &v, real4d &w, real4d &misc);

#endif


