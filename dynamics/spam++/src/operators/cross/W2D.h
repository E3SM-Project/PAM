
#ifndef _W2D_H_
#define _W2D_H_

#include "common.h"
#include "topology.h"

void YAKL_INLINE W2D_2( realArr var, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeW2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

      //x-dir
      var(0,k+ks,j+js,i+is) = -1./4. * (flux(1,k+ks,j+js  ,i+is  )
                                      + flux(1,k+ks,j+js  ,i+is-1)
                                      + flux(1,k+ks,j+js+1,i+is  )
                                      + flux(1,k+ks,j+js+1,i+is-1));

      //y-dir
      var(1,k+ks,j+js,i+is) = 1./4. * (flux(0,k+ks,j+js  ,i+is  )
                                     + flux(0,k+ks,j+js  ,i+is+1)
                                     + flux(0,k+ks,j+js-1,i+is  )
                                     + flux(0,k+ks,j+js-1,i+is+1));
  });

}

#endif
