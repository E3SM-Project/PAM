
#ifndef _DIVERGENCE_H_
#define _DIVERGENCE_H_

#include "common.h"
#include "topology.h"

template <int ndofs> void YAKL_INLINE divergence2( realArr var, realArr recon, realArr flux, Topology &topology);

#endif
