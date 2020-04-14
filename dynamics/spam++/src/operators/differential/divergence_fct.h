
#ifndef _DIVERGENCE_FCT_H_
#define _DIVERGENCE_FCT_H_

#include "common.h"
#include "topology.h"

// THIS REQUIRES SPLIT OF D into D2 K and G into L G2
// So we can do D2 PHI K and L PHI G

#define fvarx(p) flux(0,k+ks,j+js,i+is+p)*recon(l+0*ndofs,k+ks,j+js,i+is+p)
#define fvary(p) flux(1,k+ks,j+js+p,i+is)*recon(l+1*ndofs,k+ks,j+js+p,i+is)
#define fvarz(p) flux(2,k+ks+p,j+js,i+is)*recon(l+2*ndofs,k+ks+p,j+js,i+is)


template<uint ndims, uint ndofs> void YAKL_INLINE calculate_edge_fluxes2( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

      yakl::parallel_for("ComputeD2FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

        //x-dir
        var(l+0*ndofs,k+ks,j+js,i+is) = fvarx(0);

        //y-dir
        if (ndims >= 2) {
        var(l+1*ndofs,k+ks,j+js,i+is) = fvary(0);
        }

        //z-dir
        if (ndims >= 3) {
        var(l+2*ndofs,k+ks,j+js,i+is) = fvarz(0);
        }
      }
    });

}



template<uint ndims, uint ndofs> void YAKL_INLINE calculate_edge_fluxes4( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

    yakl::parallel_for("ComputeD4FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

          //x-dir
          var(l+0*ndofs,k+ks,j+js,i+is) = -1./24.*fvarx(-1) + 26./24.*fvarx(0) -1./24.*fvarx(1);

          //y-dir
          if (ndims >= 2) {
          var(l+1*ndofs,k+ks,j+js,i+is) = -1./24.*fvary(-1) + 26./24.*fvary(0) -1./24.*fvary(1);
          }

          //z-dir
          if (ndims >= 3) {
          var(l+2*ndofs,k+ks,j+js,i+is) = -1./24.*fvarz(-1) + 26./24.*fvarz(0) -1./24.*fvarz(1);
          }
    }
  });

}

template<uint ndims, uint ndofs> void YAKL_INLINE calculate_edge_fluxes6( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

      yakl::parallel_for("ComputeD6FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

            //x-dir
            var(l+0*ndofs,k+ks,j+js,i+is) = 9./1920.*fvarx(-2) - 116./1920.*fvarx(-1) + 2134./1920.*fvarx(0) - 116./1920.*fvarx(1) + 9./1920.*fvarx(2);

            //y-dir
            if (ndims >= 2) {
            var(l+1*ndofs,k+ks,j+js,i+is) = 9./1920.*fvary(-2) - 116./1920.*fvary(-1) + 2134./1920.*fvary(0) - 116./1920.*fvary(1) + 9./1920.*fvary(2);
            }

            //z-dir
            if (ndims >= 3) {
            var(l+2*ndofs,k+ks,j+js,i+is) = 9./1920.*fvarz(-2) - 116./1920.*fvarz(-1) + 2134./1920.*fvarz(0) - 116./1920.*fvarz(1) + 9./1920.*fvarz(2);
            }

      }
    });
}

//THIS IS BROKEN!!!
template<uint ndims, uint ndofs> void YAKL_INLINE calculate_edge_fluxes8( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;
    real fluxrecon0, fluxrecon1;

    yakl::parallel_for("ComputeD6FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

          //THIS IS BROKEN!!!
          //x-dir
          var(l+0*ndofs,k+ks,j+js,i+is) = 9./1920.*fvarx(-2) - 116./1920.*fvarx(-1) + 2134./1920.*fvarx(0) - 116./1920.*fvarx(1) + 9./1920.*fvarx(2);
          //THIS IS BROKEN!!!

          //y-dir
          if (ndims >= 2) {
          var(l+1*ndofs,k+ks,j+js,i+is) = 9./1920.*fvary(-2) - 116./1920.*fvary(-1) + 2134./1920.*fvary(0) - 116./1920.*fvary(1) + 9./1920.*fvary(2);
          }

          //THIS IS BROKEN!!!

          //z-dir
          if (ndims >= 3) {
          var(l+2*ndofs,k+ks,j+js,i+is) = 9./1920.*fvarz(-2) - 116./1920.*fvarz(-1) + 2134./1920.*fvarz(0) - 116./1920.*fvarz(1) + 9./1920.*fvarz(2);
          }

    }
  });
}

template<uint ndims, uint ndofs> void calculate_Mf(realArr Mf, realArr edgeflux, real dt, const Topology<ndims> &topo)
{
    int is = topo.is;
    int js = topo.js;
    int ks = topo.ks;
    real eps = 1.0e-8;

      yakl::parallel_for("ComputeMf", topo.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topo.n_cells_z, topo.n_cells_y, topo.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

            Mf(l,k+ks,j+js,i+is) = mymax(edgeflux(l+0*ndofs,k+ks,j+js,i+is+1),0.0) - mymin(edgeflux(l+0*ndofs,k+ks,j+js,i+is), 0.0);
            if (ndims >=2){
            Mf(l,k+ks,j+js,i+is) += mymax(edgeflux(l+1*ndofs,k+ks,j+js+1,i+is),0.0) - mymin(edgeflux(l+1*ndofs,k+ks,j+js,i+is), 0.0);
            }
            if (ndims ==3){
            Mf(l,k+ks,j+js,i+is) += mymax(edgeflux(l+2*ndofs,k+ks+1,j+js,i+is),0.0) - mymin(edgeflux(l+2*ndofs,k+ks,j+js,i+is), 0.0);
            }
            Mf(l,k+ks,j+js,i+is) *= dt;
            Mf(l,k+ks,j+js,i+is) += eps;
        }

});
}


template<uint ndims, uint ndofs> void calculate_phi(realArr phi, realArr q, realArr Mf, realArr edgeflux, const Topology<ndims> &topo)
{
    int is = topo.is;
    int js = topo.js;
    int ks = topo.ks;
    real upwind_param;

      yakl::parallel_for("ComputePhi", topo.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topo.n_cells_z, topo.n_cells_y, topo.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

            upwind_param = copysign(1.0, edgeflux(l+0*ndofs, k+ks, j+js, i+is));
            upwind_param = 0.5*(upwind_param + fabs(upwind_param));
            phi(l+0*ndofs,k+ks,j+js,i+is) = mymin(1., q(l,k+ks,j+js,i+is)/Mf(l,k+ks,j+js,i+is)) * (1. - upwind_param) + mymin(1., q(l,k+ks,j+js,i+is-1)/Mf(l,k+ks,j+js,i+is-1)) * upwind_param;

            if (ndims >=2) {
                upwind_param = copysign(1.0, edgeflux(l+1*ndofs, k+ks, j+js, i+is));
                upwind_param = 0.5*(upwind_param + fabs(upwind_param));
                phi(l+1*ndofs,k+ks,j+js,i+is) = mymin(1., q(l,k+ks,j+js,i+is)/Mf(l,k+ks,j+js,i+is)) * (1. - upwind_param) + mymin(1., q(l,k+ks,j+js-1,i+is)/Mf(l,k+ks,j+js-1,i+is)) * upwind_param;
            }
            if (ndims ==3) {
                upwind_param = copysign(1.0, edgeflux(l+2*ndofs, k+ks, j+js, i+is));
                upwind_param = 0.5*(upwind_param + fabs(upwind_param));
                phi(l+2*ndofs,k+ks,j+js,i+is) = mymin(1., q(l,k+ks,j+js,i+is)/Mf(l,k+ks,j+js,i+is)) * (1. - upwind_param) + mymin(1., q(l,k+ks-1,j+js,i+is)/Mf(l,k+ks-1,j+js,i+is)) * upwind_param;
        }
    }
});
}


template<uint ndims, uint ndofs> void YAKL_INLINE divergence_fct( realArr var, const realArr phi, const realArr edgeflux, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

      yakl::parallel_for("ComputeDFCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

        //x-dir
        var(l,k+ks,j+js,i+is) = phi(l+0*ndofs,k+ks,j+js,i+is+1)*edgeflux(l+0*ndofs,k+ks,j+js,i+is+1) - phi(l+0*ndofs,k+ks,j+js,i+is)*edgeflux(l+0*ndofs,k+ks,j+js,i+is);

        //y-dir
        if (ndims >= 2) {
        var(l,k+ks,j+js,i+is) += phi(l+1*ndofs,k+ks,j+js+1,i+is)*edgeflux(l+1*ndofs,k+ks,j+js+1,i+is) - phi(l+1*ndofs,k+ks,j+js,i+is)*edgeflux(l+1*ndofs,k+ks,j+js,i+is);
        }
        //z-dir
        if (ndims >= 3) {
        var(l,k+ks,j+js,i+is) += phi(l+2*ndofs,k+ks+1,j+js,i+is)*edgeflux(l+2*ndofs,k+ks+1,j+js,i+is) - phi(l+2*ndofs,k+ks,j+js,i+is)*edgeflux(l+2*ndofs,k+ks,j+js,i+is);
        }
      }
    });

}
#endif
