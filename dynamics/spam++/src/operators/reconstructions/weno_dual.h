
#ifndef _WENODUAL_H_
#define _WENODUAL_H_

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "interpolations.h"


// Indexing for storage here changes depending on the number of dimensions
// for example, in 2D 1st edge is yvar, 2nd edge is xvar
// in 3D 1st edge is zvar, 2nd edge is is yvar, 3rd edge is xvar

// star_star_sign is used to indicate that the x-dir has a minus sign for the 2D case
// IFF flux is a dual n-1 form
// otherwise it keeps the + sign

template<uint ndims, uint ndofs> void YAKL_INLINE weno1_dual_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real upwind_param;
  int star_star_sign = 1;
  if (ndims == 2) star_star_sign = -1;

  yakl::parallel_for("ComputeWENO1DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = xvar_dual(1) * (1. - upwind_param) + xvar_dual(0) * upwind_param;

      //y-dir
      if (ndims >= 2) {
      upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = yvar_dual(1) * (1. - upwind_param) + yvar_dual(0) * upwind_param;
      }
      //z-dir
      if (ndims == 3) {
      upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = zvar_dual(1) * (1. - upwind_param) + zvar_dual(0) * upwind_param;
      }
    }
  });
}

template<uint ndims, uint ndofs> void YAKL_INLINE weno3_dual_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;
  int star_star_sign = 1;
  if (ndims == 2) star_star_sign = -1;

  yakl::parallel_for("ComputeWENO3DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno3(xvar_dual(-1), xvar_dual(0), xvar_dual(1));
      var_down = interp_weno3(xvar_dual(0), xvar_dual(1), xvar_dual(2));
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno3(yvar_dual(-1), yvar_dual(0), yvar_dual(1));
        var_down = interp_weno3(yvar_dual(0), yvar_dual(1), yvar_dual(2));
        recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno3(zvar_dual(-1), zvar_dual(0), zvar_dual(1));
        var_down = interp_weno3(zvar_dual(0), zvar_dual(1), zvar_dual(2));
        recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}




template<uint ndims, uint ndofs> void YAKL_INLINE weno5_dual_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;
  int star_star_sign = 1;
  if (ndims == 2) star_star_sign = -1;

  yakl::parallel_for("ComputeWENO5DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno5(xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2));
      var_down = interp_weno5(xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3));
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno5(yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2));
        var_down = interp_weno5(yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3));
        recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno5(zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2));
        var_down = interp_weno5(zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3));
        recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}





template<uint ndims, uint ndofs> void YAKL_INLINE weno7_dual_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;
  int star_star_sign = 1;
  if (ndims == 2) star_star_sign = -1;

  yakl::parallel_for("ComputeWENO7DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno7(xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3));
      var_down = interp_weno7(xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4));
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno7(yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3));
        var_down = interp_weno7(yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4));
        recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno7(zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3));
        var_down = interp_weno7(zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4));
        recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}



template<uint ndims, uint ndofs> void YAKL_INLINE weno9_dual_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real var_up, var_down, upwind_param;
  int star_star_sign = 1;
  if (ndims == 2) star_star_sign = -1;

  yakl::parallel_for("ComputeWENO9DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
      upwind_param = 0.5*(upwind_param + fabs(upwind_param));
      var_up   = interp_weno9(xvar_dual(-4), xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4));
      var_down = interp_weno9(xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4), xvar_dual(5));
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      //y-dir
      if (ndims >= 2) {
        upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno9(yvar_dual(-4), yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4));
        var_down = interp_weno9(yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4), yvar_dual(5));
        recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

      }
      //z-dir
      if (ndims == 3) {
        upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno9(zvar_dual(-4), zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4));
        var_down = interp_weno9(zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4), zvar_dual(5));
        recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE weno11_dual_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {


    int is = topology.is;
    int js = topology.js;
    int ks = topology.ks;

    real var_up, var_down, upwind_param;
    int star_star_sign = 1;
    if (ndims == 2) star_star_sign = -1;

    yakl::parallel_for("ComputeWENO11DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
      int k, j, i;
      yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
      for (int l=0; l<ndofs; l++) {

        //x-dir
        upwind_param = copysign(1.0, flux(ndims-1, k+ks, j+js, i+is));
        upwind_param = 0.5*(upwind_param + fabs(upwind_param));
        var_up   = interp_weno11(xvar_dual(-5), xvar_dual(-4), xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4), xvar_dual(5));
        var_down = interp_weno11(xvar_dual(-4), xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4), xvar_dual(5), xvar_dual(6));
        recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

        //y-dir
        if (ndims >= 2) {
          upwind_param = copysign(1.0, star_star_sign*flux(ndims-2, k+ks, j+js, i+is));
          upwind_param = 0.5*(upwind_param + fabs(upwind_param));
          var_up   = interp_weno11(yvar_dual(-5), yvar_dual(-4), yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4), yvar_dual(5));
          var_down = interp_weno11(yvar_dual(-4), yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4), yvar_dual(5), yvar_dual(6));
          recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;

        }
        //z-dir
        if (ndims == 3) {
          upwind_param = copysign(1.0, flux(ndims-3, k+ks, j+js, i+is));
          upwind_param = 0.5*(upwind_param + fabs(upwind_param));
          var_up   = interp_weno11(zvar_dual(-5), zvar_dual(-4), zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4), zvar_dual(5));
          var_down = interp_weno11(zvar_dual(-4), zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4), zvar_dual(5), zvar_dual(6));
          recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = var_down * (1. - upwind_param) + var_up * upwind_param;
        }
      }
    });
  }


#endif
