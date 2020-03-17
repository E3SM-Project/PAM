
#ifndef _DIVERGENCE_H_
#define _DIVERGENCE_H_

#include "common.h"
#include "topology.h"

//INDEXING IS BROKEN
//ALSO SIGNS
//ORIENTATION STUFF IS TRICKY HERE WITH 2D ALSO- DUAL DIVERGENCE ACTS ON DUAL n-1 FORMS, NOT DUAL 1-FORMS ie D-grid primal winds!

template<uint ndims, uint ndofs> void YAKL_INLINE dual_divergence2( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeD2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      //x-dir
      var(l,k+ks,j+js,i+is) = flux(0,k+ks,j+js,i+is) * recon(l+0*ndofs,k+ks,j+js,i+is) - flux(0,k+ks,j+js,i+is-1) * recon(l+0*ndofs,k+ks,j+js,i+is-1);
      //y-dir
      if (ndims >= 2) {
      var(l,k+ks,j+js,i+is) += flux(1,k+ks,j+js,i+is) * recon(l+1*ndofs,k+ks,j+js+1,i+is) - flux(1,k+ks,j+js-1,i+is) * recon(l+1*ndofs,k+ks,j+js-1,i+is);
      }
      //z-dir
      if (ndims >= 3) {
      var(l,k+ks,j+js,i+is) += flux(2,k+ks,j+js,i+is) * recon(l+2*ndofs,k+ks,j+js,i+is) - flux(2,k+ks-1,j+js,i+is) * recon(l+2*ndofs,k+ks-1,j+js,i+is);
      }
    }
  });

}


template<uint ndims, uint ndofs> void YAKL_INLINE dual_divergence4( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeD4", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      //x-dir
      var(l,k+ks,j+js,i+is) = 9./8. * (flux(0,k+ks,j+js,i+is+1) * recon(l+0*ndofs,k+ks,j+js,i+is+1) - flux(0,k+ks,j+js,i+is) * recon(l+0*ndofs,k+ks,j+js,i+is))
                           - 1./24. * (flux(0,k+ks,j+js,i+is+2) * recon(l+0*ndofs,k+ks,j+js,i+is+2) - flux(0,k+ks,j+js,i+is-1) * recon(l+0*ndofs,k+ks,j+js,i+is-1));
      //y-dir
      if (ndims >= 2) {
      var(l,k+ks,j+js,i+is) += 9./8. * (flux(1,k+ks,j+js+1,i+is) * recon(l+1*ndofs,k+ks,j+js+1,i+is) - flux(1,k+ks,j+js,i+is) * recon(l+1*ndofs,k+ks,j+js,i+is))
                           - 1./24. * (flux(1,k+ks,j+js+2,i+is) * recon(l+1*ndofs,k+ks,j+js+2,i+is) - flux(1,k+ks,j+js-1,i+is) * recon(l+1*ndofs,k+ks,j+js-1,i+is));
      }
      //z-dir
      if (ndims >= 3) {
      var(l,k+ks,j+js,i+is) += 9./8. * (flux(2,k+ks+1,j+js,i+is) * recon(l+2*ndofs,k+ks+1,j+js,i+is) - flux(2,k+ks,j+js,i+is) * recon(l+2*ndofs,k+ks,j+js,i+is))
                           - 1./24. * (flux(2,k+ks+2,j+js,i+is) * recon(l+2*ndofs,k+ks+2,j+js,i+is) - flux(2,k+ks-1,j+js,i+is) * recon(l+2*ndofs,k+ks-1,j+js,i+is));
      }
    }
  });

}

template<uint ndims, uint ndofs> void YAKL_INLINE dual_divergence6( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeD6", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

    //x-dir
    var(l,k+ks,j+js,i+is) = 75./64. * (flux(0,k+ks,j+js,i+is+1) * recon(l+0*ndofs,k+ks,j+js,i+is+1) - flux(0,k+ks,j+js,i+is) * recon(l+0*ndofs,k+ks,j+js,i+is))
                         - 25./384. * (flux(0,k+ks,j+js,i+is+2) * recon(l+0*ndofs,k+ks,j+js,i+is+2) - flux(0,k+ks,j+js,i+is-1) * recon(l+0*ndofs,k+ks,j+js,i+is-1))
                         + 3./640. * (flux(0,k+ks,j+js,i+is+3) * recon(l+0*ndofs,k+ks,j+js,i+is+3) - flux(0,k+ks,j+js,i+is-2) * recon(l+0*ndofs,k+ks,j+js,i+is-2));
    //y-dir
    if (ndims >= 2) {
    var(l,k+ks,j+js,i+is) += 75./64. * (flux(1,k+ks,j+js+1,i+is) * recon(l+1*ndofs,k+ks,j+js+1,i+is) - flux(1,k+ks,j+js,i+is) * recon(l+1*ndofs,k+ks,j+js,i+is))
                         - 25./384. * (flux(1,k+ks,j+js+2,i+is) * recon(l+1*ndofs,k+ks,j+js+2,i+is) - flux(1,k+ks,j+js-1,i+is) * recon(l+1*ndofs,k+ks,j+js-1,i+is))
                         + 3./640. * (flux(1,k+ks,j+js+3,i+is) * recon(l+1*ndofs,k+ks,j+js+3,i+is) - flux(1,k+ks,j+js-2,i+is) * recon(l+1*ndofs,k+ks,j+js-2,i+is));
    }
    //z-dir
    if (ndims >= 3) {
    var(l,k+ks,j+js,i+is) += 75./64. * (flux(2,k+ks+1,j+js,i+is) * recon(l+2*ndofs,k+ks+1,j+js,i+is) - flux(2,k+ks,j+js,i+is) * recon(l+2*ndofs,k+ks,j+js,i+is))
                         - 25./384. * (flux(2,k+ks+2,j+js,i+is) * recon(l+2*ndofs,k+ks+2,j+js,i+is) - flux(2,k+ks-1,j+js,i+is) * recon(l+2*ndofs,k+ks-1,j+js,i+is))
                         + 3./640. * (flux(2,k+ks+3,j+js,i+is) * recon(l+2*ndofs,k+ks+3,j+js,i+is) - flux(2,k+ks-2,j+js,i+is) * recon(l+2*ndofs,k+ks-2,j+js,i+is));
    }
  }
});
}

template<uint ndims, uint ndofs> void YAKL_INLINE dual_divergence8( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeD8", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

    //x-dir
    var(l,k+ks,j+js,i+is) = 1225./1024. * (flux(0,k+ks,j+js,i+is+1) * recon(l+0*ndofs,k+ks,j+js,i+is+1) - flux(0,k+ks,j+js,i+is) * recon(l+0*ndofs,k+ks,j+js,i+is))
                         - 245./3072. * (flux(0,k+ks,j+js,i+is+2) * recon(l+0*ndofs,k+ks,j+js,i+is+2) - flux(0,k+ks,j+js,i+is-1) * recon(l+0*ndofs,k+ks,j+js,i+is-1))
                         + 49./5120. * (flux(0,k+ks,j+js,i+is+3) * recon(l+0*ndofs,k+ks,j+js,i+is+3) - flux(0,k+ks,j+js,i+is-2) * recon(l+0*ndofs,k+ks,j+js,i+is-2))
                         - 5./7168. * (flux(0,k+ks,j+js,i+is+4) * recon(l+0*ndofs,k+ks,j+js,i+is+4) - flux(0,k+ks,j+js,i+is-3) * recon(l+0*ndofs,k+ks,j+js,i+is-3));

    //y-dir
    if (ndims >= 2) {
    var(l,k+ks,j+js,i+is) += 1225./1024. * (flux(1,k+ks,j+js+1,i+is) * recon(l+1*ndofs,k+ks,j+js+1,i+is) - flux(1,k+ks,j+js,i+is) * recon(l+1*ndofs,k+ks,j+js,i+is))
                         - 245./3072. * (flux(1,k+ks,j+js+2,i+is) * recon(l+1*ndofs,k+ks,j+js+2,i+is) - flux(1,k+ks,j+js-1,i+is) * recon(l+1*ndofs,k+ks,j+js-1,i+is))
                         + 49./5120. * (flux(1,k+ks,j+js+3,i+is) * recon(l+1*ndofs,k+ks,j+js+3,i+is) - flux(1,k+ks,j+js-2,i+is) * recon(l+1*ndofs,k+ks,j+js-2,i+is));
                         - 5./7168. * (flux(0,k+ks,j+js+3,i+is) * recon(l+0*ndofs,k+ks,j+js+4,i+is) - flux(0,k+ks,j+js-3,i+is) * recon(l+0*ndofs,k+ks,j+js-3,i+is));
    }
    //z-dir
    if (ndims >= 3) {
    var(l,k+ks,j+js,i+is) += 1225./1024. * (flux(2,k+ks+1,j+js,i+is) * recon(l+2*ndofs,k+ks+1,j+js,i+is) - flux(2,k+ks,j+js,i+is) * recon(l+2*ndofs,k+ks,j+js,i+is))
                         - 245./3072. * (flux(2,k+ks+2,j+js,i+is) * recon(l+2*ndofs,k+ks+2,j+js,i+is) - flux(2,k+ks-1,j+js,i+is) * recon(l+2*ndofs,k+ks-1,j+js,i+is))
                         + 49./5120. * (flux(2,k+ks+3,j+js,i+is) * recon(l+2*ndofs,k+ks+3,j+js,i+is) - flux(2,k+ks-2,j+js,i+is) * recon(l+2*ndofs,k+ks-2,j+js,i+is));
                         - 5./7168. * (flux(0,k+ks+4,j+js,i+is) * recon(l+0*ndofs,k+ks+4,j+js,i+is) - flux(0,k+ks-3,j+js,i+is) * recon(l+0*ndofs,k+ks-3,j+js,i+is));
    }
  }
});


}

#endif
