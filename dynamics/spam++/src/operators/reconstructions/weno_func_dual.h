#ifndef _WENO_DUALFUNC_H_
#define _WENO_DUALFUNC_H_

#include "common.h"
#include "interpolations.h"
#include "weno_func_dual_helpers.h"


template<uint ndims, uint ndofs> void YAKL_INLINE wenodualfunc_compute_edgerecons(realArr edgerecons, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  TransformMatrices<real> trans;
  SArray<real,dual_ord,dual_tord> to_gll;
  SArray<real,dual_ord,dual_ord,dual_ord> wenoRecon;
  SArray<real,dual_hs+2> wenoIdl;
  real wenoSigma;

  // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
  trans.coefs_to_gll_lower( to_gll );
  trans.weno_sten_to_coefs(wenoRecon);
  wenoSetIdealSigma_dual(wenoIdl,wenoSigma);

  SArray<real,dual_ord> stencil;
  SArray<real,dual_tord> gllPts;

  yakl::parallel_for("ComputeWENODualFuncEdgeRecons", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      for (int ii=0; ii<dual_ord; ii++) { stencil(ii) = xvar(ii-dual_hs); }
      reconStencil_dual(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
      for (int ii=0; ii<dual_tord; ii++) { edgerecons(l*ndofs + ii + (ndims-1)*2,k+ks,j+js,i+is) = gllPts(ii); }

      if (ndims >= 2) {
        for (int ii=0; ii<dual_ord; ii++) { stencil(ii) = yvar(ii-dual_hs); }
        reconStencil_dual(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<dual_tord; ii++) { edgerecons(l*ndofs + ii + (ndims-2)*2,k+ks,j+js,i+is) = gllPts(ii); }
      }

      if (ndims == 3) {
        for (int ii=0; ii<dual_ord; ii++) { stencil(ii) = zvar(ii-dual_hs); }
        reconStencil_dual(stencil, gllPts, wenoRecon, to_gll, wenoIdl, wenoSigma);
        for (int ii=0; ii<dual_tord; ii++) { edgerecons(l*ndofs + ii + (ndims-3)*2,k+ks,j+js,i+is) = gllPts(ii); }
      }

}
});
}

template<uint ndims, uint ndofs> void YAKL_INLINE wenodualfunc_recon(realArr recon, const realArr edgerecons, const realArr flux, const Topology<ndims> &topology) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;
  int star_star_sign = 1;
  if (ndims == 2) star_star_sign = -1;

  yakl::parallel_for("ComputeWENODualFuncRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

// NOT SURE EDGERECON CHOICES HERE ARE CORRECT...
      //x-dir
      upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = edgerecons(l*ndofs + 1 + (ndims-1)*2,k+ks,j+js,i+is);
      var_down = edgerecons(l*ndofs + 0 + (ndims-1)*2,k+ks,j+js,i+is+1);
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = edgerecons(l*ndofs + 1 + (ndims-2)*2,k+ks,j+js,i+is);
        var_down = edgerecons(l*ndofs + 0 + (ndims-2)*2,k+ks,j+js+1,i+is);
        recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = edgerecons(l*ndofs + 1 + (ndims-3)*2,k+ks,j+js,i+is);
        var_down = edgerecons(l*ndofs + 0 + (ndims-3)*2,k+ks+1,j+js,i+is);
        recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }

  });
}

#endif
