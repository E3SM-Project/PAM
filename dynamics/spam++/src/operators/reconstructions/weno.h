
#ifndef _WENO_H_
#define _WENO_H_

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "interpolations.h"

template<uint ndims, uint ndofs> void YAKL_INLINE weno1_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real upwind_param;
  yakl::parallel_for("ComputeWENO1Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l+0*ndofs, k+ks, j+js, i+is) = xvar(0) * (1. - upwind_param) + xvar(-1) * upwind_param;

      //y-dir
      if (ndims >= 2) {
      upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l+1*ndofs, k+ks, j+js, i+is) = yvar(0) * (1. - upwind_param) + yvar(-1) * upwind_param;
      }
      //z-dir
      if (ndims == 3) {
      upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l+2*ndofs, k+ks, j+js, i+is) = zvar(0) * (1. - upwind_param) + zvar(-1) * upwind_param;
      }
    }
  });
}

template<uint ndims, uint ndofs> void YAKL_INLINE weno3_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;

  yakl::parallel_for("ComputeWENO3Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno3(xvar(-2), xvar(-1), xvar(0));
      var_down = interp_weno3(xvar(-1), xvar(0), xvar(1));
      recon(l+0*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno3(yvar(-2), yvar(-1), yvar(0));
        var_down = interp_weno3(yvar(-1), yvar(0), yvar(1));
        recon(l+1*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno3(zvar(-2), zvar(-1), zvar(0));
        var_down = interp_weno3(zvar(-1), zvar(0), zvar(1));
        recon(l+2*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}




template<uint ndims, uint ndofs> void YAKL_INLINE weno5_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;

  yakl::parallel_for("ComputeWENO5Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno5(xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1));
      var_down = interp_weno5(xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2));
      recon(l+0*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno5(yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1));
        var_down = interp_weno5(yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2));
        recon(l+1*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno5(zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1));
        var_down = interp_weno5(zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2));
        recon(l+2*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}





template<uint ndims, uint ndofs> void YAKL_INLINE weno7_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;

  yakl::parallel_for("ComputeWENO7Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno7(xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2));
      var_down = interp_weno7(xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3));
      recon(l+0*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno7(yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2));
        var_down = interp_weno7(yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3));
        recon(l+1*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno7(zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2));
        var_down = interp_weno7(zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3));
        recon(l+2*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}



template<uint ndims, uint ndofs> void YAKL_INLINE weno9_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;

  yakl::parallel_for("ComputeWENO9Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno9(xvar(-5), xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3));
      var_down = interp_weno9(xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3), xvar(4));
      recon(l+0*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno9(yvar(-5), yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3));
        var_down = interp_weno9(yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3), yvar(4));
        recon(l+1*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno9(zvar(-5), zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3));
        var_down = interp_weno9(zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3), zvar(4));
        recon(l+2*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE weno11_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {


    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

    real var_up, var_down, upwind_param;

    yakl::parallel_for("ComputeWENO11Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

        //x-dir
        upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno11(xvar(-6), xvar(-5), xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3), xvar(4));
        var_down = interp_weno11(xvar(-5), xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3), xvar(4), xvar(5));
        recon(l+0*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

        //y-dir
        if (ndims >= 2) {
          upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
          upwind_param = 0.5*(upwind_param + fabs(upwind_param));
          var_up   = interp_weno11(yvar(-6), yvar(-5), yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3), yvar(4));
          var_down = interp_weno11(yvar(-5), yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3), yvar(4), yvar(5));
          recon(l+1*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

        }
        //z-dir
        if (ndims == 3) {
          upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
          upwind_param = 0.5*(upwind_param + fabs(upwind_param));
          var_up   = interp_weno11(zvar(-6), zvar(-5), zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3), zvar(4));
          var_down = interp_weno11(zvar(-5), zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3), zvar(4), zvar(5));
          recon(l+2*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
        }
      }
    });
  }


#endif
