#ifndef _WENO_FUNC_H_
#define _WENO_FUNC_H_

#include "common.h"
#include "interpolations.h"
#include "weno_func_helpers.h"

template<uint ndims, uint ndofs> void YAKL_INLINE wenofunc_compute_edgerecons(realArr edgerecons, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  TransformMatrices<real> trans;
  //SArray<real,1,tord> gllWts;
  SArray<real,2,ord,tord> to_gll;
  SArray<real,3,ord,ord,ord> wenoRecon;
  SArray<real,1,hs+2> wenoIdl;
  real wenoSigma;

  // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
  trans.coefs_to_gll_lower( to_gll );
  trans.weno_sten_to_coefs(wenoRecon);
  //trans.get_gll_weights(gllWts);
  wenoSetIdealSigma(wenoIdl,wenoSigma);

  SArray<real,1,ord> stencil;
  SArray<real,1,tord> gllPts;

  yakl::parallel_for("ComputeWENOFuncEdgeRecons", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      for (int ii=0; ii<ord; ii++) { stencil(ii) = xvar(ii-hs); }
      reconStencil(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<tord; ii++) { edgerecons(l*ndofs + ii,k+ks,j+js,i+is) = gllPts(ii); }

      if (ndims >= 2) {
        for (int ii=0; ii<ord; ii++) { stencil(ii) = yvar(ii-hs); }
        reconStencil(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { edgerecons(l*ndofs + (ii+2),k+ks,j+js,i+is) = gllPts(ii); }
      }

      if (ndims == 3) {
        for (int ii=0; ii<ord; ii++) { stencil(ii) = zvar(ii-hs); }
        reconStencil(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<tord; ii++) { edgerecons(l*ndofs + (ii+4),k+ks,j+js,i+is) = gllPts(ii); }
      }

}
});
}

template<uint ndims, uint ndofs> void YAKL_INLINE wenofunc_recon(realArr recon, const realArr edgerecons, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;

  yakl::parallel_for("ComputeWENOFuncRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = edgerecons(l*ndofs + 1,k+ks,j+js,i+is-1);
      var_down = edgerecons(l*ndofs + 0,k+ks,j+js,i+is);
      recon(l+0*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = edgerecons(l*ndofs + 1+2,k+ks,j+js-1,i+is);
        var_down = edgerecons(l*ndofs + 0+2,k+ks,j+js,i+is);
        recon(l+1*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = edgerecons(l*ndofs + 1+4,k+ks-1,j+js,i+is);
        var_down = edgerecons(l*ndofs + 0+4,k+ks,j+js,i+is);
        recon(l+2*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}

#endif
