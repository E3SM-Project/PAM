

#include "divergence.h"


//REALLY NDOFS HERE SHOULD BE A COMPILE TIME CONSTANT...

template<int ndofs, int ndims> void YAKL_INLINE divergence2( realArr var, realArr recon, realArr flux, Topology &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;


  // PARALLEL LOOP BOUNDS?
    yakl::parallel_for("ComputeD2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      //x-dir
      var(l,k+ks,j+js,i+is) = flux(0,k+ks,j+js,i+is+1) * recon(l+0*ndofs,k+ks,j+js,i+is+1) - flux(0,k+ks,j+js,i+is) * recon(l+0*ndofs,k+ks,j+js,i+is);
      //y-dir
      if (ndims >= 2) {
      var(l,k+ks,j+js,i+is) = flux(1,k+ks,j+js+1,i+is) * recon(l+1*ndofs,k+ks,j+js+1,i+is) - flux(1,k+ks,j+js,i+is) * recon(l+1*ndofs,k+ks,j+js,i+is);
      }
      //z-dir
      if (ndims >= 3) {
      var(l,k+ks,j+js,i+is) = flux(2,k+ks+1,j+js,i+is) * recon(l+2*ndofs,k+ks+1,j+js,i+is) - flux(2,k+ks,j+js,i+is) * recon(l+2*ndofs,k+ks,j+js,i+is);
      }
    }
  });

}
}
