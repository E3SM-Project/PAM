
#ifndef _Q2D_H_
#define _Q2D_H_

#include "common.h"
#include "topology.h"

// Q = 1/2 * (q W + W q)
void YAKL_INLINE Q2D_2_add( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeQ2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

    //x-dir
    var(0,k+ks,j+js,i+is) += -1./4. * ((recon(0,k+ks,j+js,i+is) + recon(1,k+ks,j+js  ,i+is  ))/2. * flux(1,k+ks,j+js,  i+is  )
                                    + (recon(0,k+ks,j+js,i+is) + recon(1,k+ks,j+js  ,i+is-1))/2. * flux(1,k+ks,j+js,  i+is-1)
                                    + (recon(0,k+ks,j+js,i+is) + recon(1,k+ks,j+js+1,i+is  ))/2. * flux(1,k+ks,j+js+1,i+is  )
                                    + (recon(0,k+ks,j+js,i+is) + recon(1,k+ks,j+js+1,i+is-1))/2. * flux(1,k+ks,j+js+1,i+is-1));

    //y-dir
    var(1,k+ks,j+js,i+is) +=  1./4. * ((recon(1,k+ks,j+js,i+is) + recon(0,k+ks,j+js  ,i+is  ))/2. * flux(0,k+ks,j+js  ,i+is  )
                                   +  (recon(1,k+ks,j+js,i+is) + recon(0,k+ks,j+js  ,i+is+1))/2. * flux(0,k+ks,j+js  ,i+is+1)
                                   +  (recon(1,k+ks,j+js,i+is) + recon(0,k+ks,j+js-1,i+is  ))/2. * flux(0,k+ks,j+js-1,i+is  )
                                   +  (recon(1,k+ks,j+js,i+is) + recon(0,k+ks,j+js-1,i+is+1))/2. * flux(0,k+ks,j+js-1,i+is+1));
});

}


// Q = q W
void YAKL_INLINE Q2D_nonEC_2_add( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeQ2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);

   //x-dir
   var(0,k+ks,j+js,i+is) += -1./4. * recon(0,k+ks,j+js,i+is) * (flux(1,k+ks,j+js  ,i+is  )
                                                            + flux(1,k+ks,j+js  ,i+is-1)
                                                            + flux(1,k+ks,j+js+1,i+is  )
                                                            + flux(1,k+ks,j+js+1,i+is-1));

   //y-dir
   var(1,k+ks,j+js,i+is) += 1./4. * recon(1,k+ks,j+js,i+is) * (flux(0,k+ks,j+js  ,i+is  )
                                                            + flux(0,k+ks,j+js  ,i+is+1)
                                                            + flux(0,k+ks,j+js-1,i+is  )
                                                            + flux(0,k+ks,j+js-1,i+is+1));

});

}

#endif
