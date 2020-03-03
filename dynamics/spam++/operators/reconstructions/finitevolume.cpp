

#include "finitevolume.h"

//REALLY NDOFS HERE SHOULD BE A COMPILE TIME CONSTANT...

void YAKL_INLINE fv1_recon( realArr recon, realArr var, int ndofs, Topology &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  // FIX THIS TO REFER TO PRIMAL CELLS, ETC.

  // PARALLEL LOOP BOUNDS?
  yakl::parallel_for("ComputeFV1Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*nqdofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1))/2.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*nqdofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is))/2.0;
      }
      //z-dir
      if (ndims >= 3) {
      recon(l+2*nqdofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is))/2.0;
      }
    }
  });
}
