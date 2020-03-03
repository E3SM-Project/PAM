
#ifndef _FINITEVOLUME_H_
#define _FINITEVOLUME_H_

#include "common.h"
#include "topology.h"

void YAKL_INLINE fv1_recon( realArr recon, realArr var, int ndofs, Topology &topology);

#endif
