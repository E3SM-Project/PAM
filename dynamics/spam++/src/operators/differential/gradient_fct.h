
#ifndef _DIVERGENCE_FCT_H_
#define _DIVERGENCE_FCT_H_

#include "common.h"
#include "topology.h"

// THIS REQUIRES SPLIT OF D into D2 K and G into L G2
// So we can do D2 PHI K and L PHI G

#define fvarx(p) flux(0,k+ks,j+js,i+is+p)*recon(l+0*ndofs,k+ks,j+js,i+is+p)
#define fvary(p) flux(1,k+ks,j+js+p,i+is)*recon(l+1*ndofs,k+ks,j+js+p,i+is)
#define fvarz(p) flux(2,k+ks+p,j+js,i+is)*recon(l+2*ndofs,k+ks+p,j+js,i+is)


template<uint ndims, uint ndofs> void YAKL_INLINE divergence_fct_2( realArr var, const realArr recon, const realArr flux, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;
    real fluxrecon0, fluxrecon1;

      yakl::parallel_for("ComputeD2FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

        //x-dir
        fluxrecon1 = fvarx(1);
        fluxrecon0 = fvarx(0);
        var(l,k+ks,j+js,i+is) = phi(l+0*ndofs,k+ks,j+js,i+is+1)*fluxrecon1 - phi(l+0*ndofs,k+ks,j+js,i+is)*fluxrecon0;

        //y-dir
        if (ndims >= 2) {
        fluxrecon1 = fvary(1);
        fluxrecon0 = fvary(0);
        var(l,k+ks,j+js,i+is) += phi(l+1*ndofs,k+ks,j+js+1,i+is)*fluxrecon1 - phi(l+1*ndofs,k+ks,j+js,i+is)*fluxrecon0;
        }
        //z-dir
        if (ndims >= 3) {
        fluxrecon1 = fvarz(1);
        fluxrecon0 = fvarz(0);
        var(l,k+ks,j+js,i+is) += phi(l+2*ndofs,k+ks+1,j+js,i+is)*fluxrecon1 - phi(l+2*ndofs,k+ks,j+js,i+is)*fluxrecon0;
        }
      }
    });

}



template<uint ndims, uint ndofs> void YAKL_INLINE divergence_fct_4( realArr var, const realArr recon, const realArr flux, const realArr phi, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;
  real fluxrecon0, fluxrecon1;

    yakl::parallel_for("ComputeD4FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

      //x-dir
      fluxrecon1 = -1./24.*fvarx(0) + 26./24.*fvarx(1) -1./24.*fvarx(2);
      fluxrecon0 = -1./24.*fvarx(-1) + 26./24.*fvarx(0) -1./24.*fvarx(1);
      var(l,k+ks,j+js,i+is) = phi(l+0*ndofs,k+ks,j+js,i+is+1)*fluxrecon1 - phi(l+0*ndofs,k+ks,j+js,i+is)*fluxrecon0;

      //y-dir
      if (ndims >= 2) {
      fluxrecon1 = -1./24.*fvary(0) + 26./24.*fvary(1) -1./24.*fvary(2);
      fluxrecon0 = -1./24.*fvary(-1) + 26./24.*fvary(0) -1./24.*fvary(1);
      var(l,k+ks,j+js,i+is) += phi(l+1*ndofs,k+ks,j+js+1,i+is)*fluxrecon1 - phi(l+1*ndofs,k+ks,j+js,i+is)*fluxrecon0;
      }
      //z-dir
      if (ndims >= 3) {
      fluxrecon1 = -1./24.*fvarz(0) + 26./24.*fvarz(1) -1./24.*fvarz(2);
      fluxrecon0 = -1./24.*fvarz(-1) + 26./24.*fvarz(0) -1./24.*fvarz(1);
      var(l,k+ks,j+js,i+is) += phi(l+2*ndofs,k+ks+1,j+js,i+is)*fluxrecon1 - phi(l+2*ndofs,k+ks,j+js,i+is)*fluxrecon0;
      }
    }
  });

}

template<uint ndims, uint ndofs> void YAKL_INLINE divergence_fct_6( realArr var, const realArr recon, const realArr flux, const realArr phi, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;
    real fluxrecon0, fluxrecon1;

      yakl::parallel_for("ComputeD6FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

        //x-dir
        fluxrecon1 = 9./1920.*fvarx(-1) - 116./1920.*fvarx(0) + 2134./1920.*fvarx(1) - 116./1920.*fvarx(2) + 9./1920.*fvarx(3);
        fluxrecon0 = 9./1920.*fvarx(-2) - 116./1920.*fvarx(-1) + 2134./1920.*fvarx(0) - 116./1920.*fvarx(1) + 9./1920.*fvarx(2);
        var(l,k+ks,j+js,i+is) = phi(l+0*ndofs,k+ks,j+js,i+is+1)*fluxrecon1 - phi(l+0*ndofs,k+ks,j+js,i+is)*fluxrecon0;

        //y-dir
        if (ndims >= 2) {
        fluxrecon1 = 9./1920.*fvary(-1) - 116./1920.*fvary(0) + 2134./1920.*fvary(1) - 116./1920.*fvary(2) + 9./1920.*fvary(3);
        fluxrecon0 = 9./1920.*fvary(-2) - 116./1920.*fvary(-1) + 2134./1920.*fvary(0) - 116./1920.*fvary(1) + 9./1920.*fvary(2);
        var(l,k+ks,j+js,i+is) += phi(l+1*ndofs,k+ks,j+js+1,i+is)*fluxrecon1 - phi(l+1*ndofs,k+ks,j+js,i+is)*fluxrecon0;
        }
        //z-dir
        if (ndims >= 3) {
        fluxrecon1 = 9./1920.*fvarz(-1) - 116./1920.*fvarz(0) + 2134./1920.*fvarz(1) - 116./1920.*fvarz(2) + 9./1920.*fvarz(3);
        fluxrecon0 = 9./1920.*fvarz(-2) - 116./1920.*fvarz(-1) + 2134./1920.*fvarz(0) - 116./1920.*fvarz(1) + 9./1920.*fvarz(2);
        var(l,k+ks,j+js,i+is) += phi(l+2*ndofs,k+ks+1,j+js,i+is)*fluxrecon1 - phi(l+2*ndofs,k+ks,j+js,i+is)*fluxrecon0;
        }
      }
    });
}

//THIS IS BROKEN!!!
template<uint ndims, uint ndofs> void YAKL_INLINE divergence_fct_8( realArr var, const realArr recon, const realArr flux, const realArr phi, const Topology<ndims> &topology) {

    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;
    real fluxrecon0, fluxrecon1;

      yakl::parallel_for("ComputeD8FCT", topology.n_cells, YAKL_LAMBDA (int iGlob) {
        int k, j, i;
        yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
        for (int l=0; l<ndofs; l++) {

//THIS IS BROKEN!!!
        //x-dir
        fluxrecon1 = 9./1920.*fvarx(-1) - 116./1920.*fvarx(0) + 2134./1920.*fvarx(1) - 116./1920*fvarx(2) + 9./1920.*fvarx(3);
        fluxrecon0 = 9./1920.*fvarx(-2) - 116./1920.*fvarx(-1) + 2134./1920.*fvarx(0) - 116./1920*fvarx(1) + 9./1920.*fvarx(2);
        var(l,k+ks,j+js,i+is) = phi(l+0*ndofs,k+ks,j+js,i+is+1)*fluxrecon1 - phi(l+0*ndofs,k+ks,j+js,i+is)*fluxrecon0;

//THIS IS BROKEN!!!
        //y-dir
        if (ndims >= 2) {
        fluxrecon1 = 9./1920.*fvary(-1) - 116./1920.*fvary(0) + 2134./1920.*fvary(1) - 116./1920*fvary(2) + 9./1920.*fvary(3);
        fluxrecon0 = 9./1920.*fvary(-2) - 116./1920.*fvary(-1) + 2134./1920.*fvary(0) - 116./1920*fvary(1) + 9./1920.*fvary(2);
        var(l,k+ks,j+js,i+is) += phi(l+1*ndofs,k+ks,j+js+1,i+is)*fluxrecon1 - phi(l+1*ndofs,k+ks,j+js,i+is)*fluxrecon0;
        }

//THIS IS BROKEN!!!
        //z-dir
        if (ndims >= 3) {
        fluxrecon1 = 9./1920.*fvarz(-1) - 116./1920.*fvarz(0) + 2134./1920.*fvarz(1) - 116./1920*fvarz(2) + 9./1920.*fvarz(3);
        fluxrecon0 = 9./1920.*fvarz(-2) - 116./1920.*fvarz(-1) + 2134./1920.*fvarz(0) - 116./1920*fvarz(1) + 9./1920.*fvarz(2);
        var(l,k+ks,j+js,i+is) += phi(l+2*ndofs,k+ks+1,j+js,i+is)*fluxrecon1 - phi(l+2*ndofs,k+ks,j+js,i+is)*fluxrecon0;
        }
      }
    });
}
#endif
