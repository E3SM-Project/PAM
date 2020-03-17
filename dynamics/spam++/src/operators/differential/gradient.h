
#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include "common.h"
#include "topology.h"

template<uint ndims, uint ndofs> void YAKL_INLINE gradient2( realArr var, const realArr recon, const realArr dens, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeG2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      //x-dir
      var(l+0*ndofs,k+ks,j+js,i+is) = recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js,i+is-1));
      //y-dir
      if (ndims >= 2) {
      var(l+1*ndofs,k+ks,j+js,i+is) = recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js-1,i+is);
      }
      //z-dir
      if (ndims >= 3) {
      var(l+2*ndofs,k+ks,j+js,i+is) = recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks-1,j+js,i+is));
      }
    }
  });

}


template<uint ndims, uint ndofs> void YAKL_INLINE gradient4( realArr var, const realArr recon, const realArr dens, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeG4", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      //x-dir
      var(l+0*ndofs,k+ks,j+js,i+is) = 9./8. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js,i+is-1))
                                   - 1./24. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is+1) - dens(l,k+ks,j+js,i+is-2));
      //y-dir
      if (ndims >= 2) {
      var(l+1*ndofs,k+ks,j+js,i+is) = 9./8. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js-1,i+is))
                                   - 1./24. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js+1,i+is) - dens(l,k+ks,j+js-2,i+is));
      }
      //z-dir
      if (ndims >= 3) {
      var(l+2*ndofs,k+ks,j+js,i+is) = 9./8. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks-1,j+js,i+is  ))
                                   - 1./24. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks+1,j+js,i+is) - dens(l,k+ks-2,j+js,i+is));
      }
    }
  });

}

template<uint ndims, uint ndofs> void YAKL_INLINE gradient6( realArr var, const realArr recon, const realArr dens, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeG6", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

    //x-dir
    var(l+0*ndofs,k+ks,j+js,i+is) = 75./64. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js,i+is-1))
                                 - 25./384. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is+1) - dens(l,k+ks,j+js,i+is-2))
                                 + 3./640. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is+2) - dens(l,k+ks,j+js,i+is-3));

    //y-dir
    if (ndims >= 2) {
    var(l+1*ndofs,k+ks,j+js,i+is) = 75./64. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js-1,i+is))
                                 - 25./384. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js+1,i+is) - dens(l,k+ks,j+js-2,i+is))
                                 + 3./640. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js+2,i+is) - dens(l,k+ks,j+js-3,i+is));
    }
    //z-dir
    if (ndims >= 3) {
    var(l+2*ndofs,k+ks,j+js,i+is) = 75./64. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks-1,j+js,i+is  ))
                                 - 25./384. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks+1,j+js,i+is) - dens(l,k+ks-2,j+js,i+is))
                                 + 3./640. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks+2,j+js,i+is) - dens(l,k+ks-3,j+js,i+is));
    }
  }
});

}

template<uint ndims, uint ndofs> void YAKL_INLINE gradient8( realArr var, const realArr recon, const realArr dens, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeG8", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

    //x-dir
    var(l+0*ndofs,k+ks,j+js,i+is) = 1225./1024. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js,i+is-1))
                                 - 245./3072. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is+1) - dens(l,k+ks,j+js,i+is-2))
                                 + 49./5120. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is+2) - dens(l,k+ks,j+js,i+is-3))
                                 - 5./7168. * recon(l+0*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is+3) - dens(l,k+ks,j+js,i+is-4));

    //y-dir
    if (ndims >= 2) {
    var(l+1*ndofs,k+ks,j+js,i+is) = 1225./1024. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks,j+js-1,i+is))
                                 - 245./3072. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js+1,i+is) - dens(l,k+ks,j+js-2,i+is))
                                 + 49./5120. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js+2,i+is) - dens(l,k+ks,j+js-3,i+is))
                                 - 5./7168. * recon(l+1*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js+3,i+is) - dens(l,k+ks,j+js-4,i+is));
    }
    //z-dir
    if (ndims >= 3) {
    var(l+2*ndofs,k+ks,j+js,i+is) = 1225./1024. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks,j+js,i+is) - dens(l,k+ks-1,j+js,i+is  ))
                                 - 245./3072. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks+1,j+js,i+is) - dens(l,k+ks-2,j+js,i+is))
                                 + 49./5120. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks+2,j+js,i+is) - dens(l,k+ks-3,j+js,i+is))
                                 - 5./7168. * recon(l+2*ndofs,k+ks,j+js,i+is) * (dens(l,k+ks+3,j+js,i+is) - dens(l,k+ks-4,j+js,i+is));
    }
  }
});

}

#endif
