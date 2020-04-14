
#ifndef _CURL_H_
#define _CURL_H_

#include "common.h"
#include "topology.h"

template<uint ndofs> void YAKL_INLINE curl2D_2( realArr var, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeC2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      var(l,k+ks,j+js,i+is) = (flux(l+1*ndofs,k+ks,j+js,i+is) - flux(l+0*ndofs,k+ks,j+js,i+is) - flux(l+1*ndofs,k+ks,j+js,i+is-1) + flux(l+0*ndofs,k+ks,j+js-1,i+is));
  }
  });

}

template<uint ndofs> void YAKL_INLINE weighted_curl2D_2( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeC2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      var(l,k+ks,j+js,i+is) = (recon(l+1*ndofs,k+ks,j+js,i+is)*flux(l+1*ndofs,k+ks,j+js,i+is)
                             - recon(l+0*ndofs,k+ks,j+js,i+is)*flux(l+0*ndofs,k+ks,j+js,i+is)
                             - recon(l+1*ndofs,k+ks,j+js,i+is-1)*flux(l+1*ndofs,k+ks,j+js,i+is-1)
                             + recon(l+0*ndofs,k+ks,j+js-1,i+is)*flux(l+0*ndofs,k+ks,j+js-1,i+is));
}
  });

}
#endif
