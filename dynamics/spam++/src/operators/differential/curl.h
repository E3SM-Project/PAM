
#ifndef _CURL_H_
#define _CURL_H_

#include "common.h"
#include "topology.h"

void YAKL_INLINE curl2D_2( realArr var, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeC2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

      var(0,k+ks,j+js,i+is) = (flux(1,k+ks,j+js,i+is) - flux(0,k+ks,j+js,i+is) - flux(1,k+ks,j+js,i+is-1) + flux(0,k+ks,j+js-1,i+is));

  });

}

#endif
