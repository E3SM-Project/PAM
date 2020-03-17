
#ifndef _CFVDUAL_H_
#define _CFVDUAL_H_

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "interpolations.h"

// Indexing for storage here changes depending on the number of dimensions
// for example, in 2D 1st edge is yvar, 2nd edge is xvar
// in 3D 1st edge is zvar, 2nd edge is is yvar, 3rd edge is xvar

template<uint ndims, uint ndofs> void YAKL_INLINE cfv2_dual_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV2DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = interp_2(xvar_dual(0), xvar_dual(1));
      //y-dir
      if (ndims >= 2) {
      recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = interp_2(yvar_dual(0), yvar_dual(1));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = interp_2(zvar_dual(0), zvar_dual(1));
      }
    }
  });
}

template<uint ndims, uint ndofs> void YAKL_INLINE cfv4_dual_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV4DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = interp_4(xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2));
      //y-dir
      if (ndims >= 2) {
      recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = interp_4(yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = interp_4(zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2));
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE cfv6_dual_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV6DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = interp_6(xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3));

      //y-dir
      if (ndims >= 2) {
      recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = interp_6(yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = interp_6(zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3));
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE cfv8_dual_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV8DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = interp_8(xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4));
      //y-dir
      if (ndims >= 2) {
      recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = interp_8(yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = interp_8(zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4));
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE cfv10_dual_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV10DualRecon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+(ndims-1)*ndofs, k+ks, j+js, i+is) = interp_10(xvar_dual(-4), xvar_dual(-3), xvar_dual(-2), xvar_dual(-1), xvar_dual(0), xvar_dual(1), xvar_dual(2), xvar_dual(3), xvar_dual(4), xvar_dual(5));
      //y-dir
      if (ndims >= 2) {
      recon(l+(ndims-2)*ndofs, k+ks, j+js, i+is) = interp_10(yvar_dual(-4), yvar_dual(-3), yvar_dual(-2), yvar_dual(-1), yvar_dual(0), yvar_dual(1), yvar_dual(2), yvar_dual(3), yvar_dual(4), yvar_dual(5));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+(ndims-3)*ndofs, k+ks, j+js, i+is) = interp_10(zvar_dual(-4), zvar_dual(-3), zvar_dual(-2), zvar_dual(-1), zvar_dual(0), zvar_dual(1), zvar_dual(2), zvar_dual(3), zvar_dual(4), zvar_dual(5));
      }
    }
  });
}

#endif
