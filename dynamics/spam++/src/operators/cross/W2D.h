
#ifndef _W2D_H_
#define _W2D_H_

#include "common.h"
#include "topology.h"

template<uint ndofs> void YAKL_INLINE W2D_2( realArr var, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeW2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

      for (int l=0; l<ndofs; l++) {
      //x-dir
      var(l+0*ndofs,k+ks,j+js,i+is) = -1./4. * (flux(l+1*ndofs,k+ks,j+js  ,i+is  )
                                      + flux(l+1*ndofs,k+ks,j+js  ,i+is-1)
                                      + flux(l+1*ndofs,k+ks,j+js+1,i+is  )
                                      + flux(l+1*ndofs,k+ks,j+js+1,i+is-1));

      //y-dir
      var(l+1*ndofs,k+ks,j+js,i+is) = 1./4. * (flux(l+0*ndofs,k+ks,j+js  ,i+is  )
                                     + flux(l+0*ndofs,k+ks,j+js  ,i+is+1)
                                     + flux(l+0*ndofs,k+ks,j+js-1,i+is  )
                                     + flux(l+0*ndofs,k+ks,j+js-1,i+is+1));
                                 }
  });

}

#endif
