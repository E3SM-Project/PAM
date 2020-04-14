
#ifndef _Q2D_H_
#define _Q2D_H_

#include "common.h"
#include "topology.h"

// Q = 1/2 * (q W + W q)
template<uint ndofs> void YAKL_INLINE Q2D_2( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeQ2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

    //x-dir
    var(l+0*ndofs,k+ks,j+js,i+is) = -1./4. * ((recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js  ,i+is  ))/2. * flux(l+1*ndofs,k+ks,j+js,  i+is  )
                                    + (recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js  ,i+is-1))/2. * flux(l+1*ndofs,k+ks,j+js,  i+is-1)
                                    + (recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js+1,i+is  ))/2. * flux(l+1*ndofs,k+ks,j+js+1,i+is  )
                                    + (recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js+1,i+is-1))/2. * flux(l+1*ndofs,k+ks,j+js+1,i+is-1));

    //y-dir
    var(l+1*ndofs,k+ks,j+js,i+is) =  1./4. * ((recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js  ,i+is  ))/2. * flux(l+0*ndofs,k+ks,j+js  ,i+is  )
                                   +  (recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js  ,i+is+1))/2. * flux(l+0*ndofs,k+ks,j+js  ,i+is+1)
                                   +  (recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js-1,i+is  ))/2. * flux(l+0*ndofs,k+ks,j+js-1,i+is  )
                                   +  (recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js-1,i+is+1))/2. * flux(l+0*ndofs,k+ks,j+js-1,i+is+1));
}
});

}

// Q = 1/2 * (q W + W q)
template<uint ndofs> void YAKL_INLINE Q2D_2_add( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeQ2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

    //x-dir
    var(l+0*ndofs,k+ks,j+js,i+is) += -1./4. * ((recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js  ,i+is  ))/2. * flux(l+1*ndofs,k+ks,j+js,  i+is  )
                                    + (recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js  ,i+is-1))/2. * flux(l+1*ndofs,k+ks,j+js,  i+is-1)
                                    + (recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js+1,i+is  ))/2. * flux(l+1*ndofs,k+ks,j+js+1,i+is  )
                                    + (recon(l+0*ndofs,k+ks,j+js,i+is) + recon(l+1*ndofs,k+ks,j+js+1,i+is-1))/2. * flux(l+1*ndofs,k+ks,j+js+1,i+is-1));

    //y-dir
    var(l+1*ndofs,k+ks,j+js,i+is) +=  1./4. * ((recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js  ,i+is  ))/2. * flux(l+0*ndofs,k+ks,j+js  ,i+is  )
                                   +  (recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js  ,i+is+1))/2. * flux(l+0*ndofs,k+ks,j+js  ,i+is+1)
                                   +  (recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js-1,i+is  ))/2. * flux(l+0*ndofs,k+ks,j+js-1,i+is  )
                                   +  (recon(l+1*ndofs,k+ks,j+js,i+is) + recon(l+0*ndofs,k+ks,j+js-1,i+is+1))/2. * flux(l+0*ndofs,k+ks,j+js-1,i+is+1));
}
});

}

// Q = q W
template<uint ndofs> void YAKL_INLINE Q2D_nonEC_2( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeQ2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

   //x-dir
   var(0,k+ks,j+js,i+is) = -1./4. * recon(l+0*ndofs,k+ks,j+js,i+is) * (flux(l+1*ndofs,k+ks,j+js  ,i+is  )
                                                            + flux(l+1*ndofs,k+ks,j+js  ,i+is-1)
                                                            + flux(l+1*ndofs,k+ks,j+js+1,i+is  )
                                                            + flux(l+1*ndofs,k+ks,j+js+1,i+is-1));

   //y-dir
   var(l+1*ndofs,k+ks,j+js,i+is) = 1./4. * recon(l+1*ndofs,k+ks,j+js,i+is) * (flux(l+0*ndofs,k+ks,j+js  ,i+is  )
                                                            + flux(l+0*ndofs,k+ks,j+js  ,i+is+1)
                                                            + flux(l+0*ndofs,k+ks,j+js-1,i+is  )
                                                            + flux(l+0*ndofs,k+ks,j+js-1,i+is+1));
}
});

}

// Q = q W
template<uint ndofs> void YAKL_INLINE Q2D_nonEC_2_add( realArr var, const realArr recon, const realArr flux, const Topology<2> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeQ2D_2", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

   //x-dir
   var(l+0*ndofs,k+ks,j+js,i+is) += -1./4. * recon(l+0*ndofs,k+ks,j+js,i+is) * (flux(l+1*ndofs,k+ks,j+js  ,i+is  )
                                                            + flux(l+1*ndofs,k+ks,j+js  ,i+is-1)
                                                            + flux(l+1*ndofs,k+ks,j+js+1,i+is  )
                                                            + flux(l+1*ndofs,k+ks,j+js+1,i+is-1));

   //y-dir
   var(l+1*ndofs,k+ks,j+js,i+is) += 1./4. * recon(l+1*ndofs,k+ks,j+js,i+is) * (flux(l+0*ndofs,k+ks,j+js  ,i+is  )
                                                            + flux(l+0*ndofs,k+ks,j+js  ,i+is+1)
                                                            + flux(l+0*ndofs,k+ks,j+js-1,i+is  )
                                                            + flux(l+0*ndofs,k+ks,j+js-1,i+is+1));
}
});

}


template<uint ndofs> void calculate_Mf_dual(realArr Mf, realArr edgeflux, real dt, const Topology<2> &topo)
{
    int is = topo.is;
    int js = topo.js;
    int ks = topo.ks;
    real eps = 1.0e-8;

      yakl::parallel_for("ComputeMfDual", topo.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topo.n_cells_z, topo.n_cells_y, topo.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

            Mf(l,k+ks,j+js,i+is) = mymax(edgeflux(l+0*ndofs,k+ks,j+js-1,i+is),0.0) - mymin(edgeflux(l+0*ndofs,k+ks,j+js,i+is), 0.0);
            Mf(l,k+ks,j+js,i+is) += mymax(edgeflux(l+1*ndofs,k+ks,j+js,i+is),0.0) - mymin(edgeflux(l+1*ndofs,k+ks,j+js,i+is-1), 0.0);
            Mf(l,k+ks,j+js,i+is) *= dt;
            Mf(l,k+ks,j+js,i+is) += eps;
        }

});
}

//FIX- THIS IS SLIGHTLY BROKEN I THINK...
template<uint ndofs> void calculate_phi_dual(realArr phi, realArr q, realArr Mf, realArr edgeflux, const Topology<2> &topo)
{
    int is = topo.is;
    int js = topo.js;
    int ks = topo.ks;
    real upwind_param;

      yakl::parallel_for("ComputePhiDual", topo.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topo.n_cells_z, topo.n_cells_y, topo.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

            upwind_param = copysign(1.0, edgeflux(l+0*ndofs, k+ks, j+js, i+is));
            upwind_param = 0.5*(upwind_param + fabs(upwind_param));
            phi(l+0*ndofs,k+ks,j+js,i+is) = mymin(1., q(l,k+ks,j+js,i+is)/Mf(l,k+ks,j+js,i+is)) * (1. - upwind_param) + mymin(1., q(l,k+ks,j+js,i+is-1)/Mf(l,k+ks,j+js,i+is-1)) * upwind_param;

            upwind_param = copysign(1.0, -edgeflux(l+1*ndofs, k+ks, j+js, i+is));
            upwind_param = 0.5*(upwind_param + fabs(upwind_param));
            phi(l+1*ndofs,k+ks,j+js,i+is) = mymin(1., q(l,k+ks,j+js,i+is)/Mf(l,k+ks,j+js,i+is)) * (1. - upwind_param) + mymin(1., q(l,k+ks,j+js-1,i+is)/Mf(l,k+ks,j+js-1,i+is)) * upwind_param;

    }
});
}
#endif
