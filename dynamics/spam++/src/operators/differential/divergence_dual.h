
#ifndef _DIVERGENCE_H_
#define _DIVERGENCE_H_

#include "common.h"
#include "topology.h"

// IN 2D this is just negative of Curl operator...

template<uint ndims, uint ndofs> void YAKL_INLINE dual_divergence2( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeD2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      if (ndims == 1) {
          var(l,k+ks,j+js,i+is) = flux(ndims-1,k+ks,j+js,i+is) * recon(l+(ndims-1)*ndofs,k+ks,j+js,i+is) - flux(ndims-1,k+ks,j+js,i+is-1) * recon(l+(ndims-1)*ndofs,k+ks,j+js,i+is-1);
      }

      if (ndims == 2) {
      var(l,k+ks,j+js,i+is) = - flux(ndims-1,k+ks,j+js,i+is) * recon(l+(ndims-1)*ndofs,k+ks,j+js,i+is) + flux(ndims-1,k+ks,j+js,i+is-1) * recon(l+(ndims-1)*ndofs,k+ks,j+js,i+is-1)
                              + flux(ndims-2,k+ks,j+js,i+is) * recon(l+(ndims-2)*ndofs,k+ks,j+js,i+is) - flux(ndims-2,k+ks,j+js-1,i+is) * recon(l+(ndims-2)*ndofs,k+ks,j+js-1,i+is);
      }

      if (ndims == 3) {
          var(l,k+ks,j+js,i+is) =   flux(ndims-1,k+ks,j+js,i+is) * recon(l+(ndims-1)*ndofs,k+ks,j+js,i+is) - flux(ndims-1,k+ks,j+js,i+is-1) * recon(l+(ndims-1)*ndofs,k+ks,j+js,i+is-1)
                                  + flux(ndims-2,k+ks,j+js,i+is) * recon(l+(ndims-2)*ndofs,k+ks,j+js,i+is) - flux(ndims-2,k+ks,j+js-1,i+is) * recon(l+(ndims-2)*ndofs,k+ks,j+js-1,i+is);
                                  + flux(ndims-3,k+ks,j+js,i+is) * recon(l+(ndims-3)*ndofs,k+ks,j+js,i+is) - flux(ndims-3,k+ks-1,j+js,i+is) * recon(l+(ndims-3)*ndofs,k+ks-1,j+js,i+is);
      }
    }
  });

}

#endif
