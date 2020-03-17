
#ifndef _CFV_H_
#define _CFV_H_

#include "common.h"
#include "topology.h"
#include "geometry.h"
#include "interpolations.h"

template<uint ndims, uint ndofs> void YAKL_INLINE cfv2_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV2Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = interp_2(xvar(-1), xvar(0));
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = interp_2(yvar(-1), yvar(0));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = interp_2(zvar(-1), zvar(0));
      }
    }
  });
}

template<uint ndims, uint ndofs> void YAKL_INLINE cfv4_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV4Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = interp_4(xvar(-2), xvar(-1), xvar(0), xvar(1));
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = interp_4(yvar(-2), yvar(-1), yvar(0), yvar(1));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = interp_4(zvar(-2), zvar(-1), zvar(0), zvar(1));
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE cfv6_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV6Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = interp_6(xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2));

      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = interp_6(yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = interp_6(zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2));
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE cfv8_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV8Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = interp_8(xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3));
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = interp_8(yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = interp_8(zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3));
      }
    }
  });
}


template<uint ndims, uint ndofs> void YAKL_INLINE cfv10_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV10Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = interp_10(xvar(-5), xvar(-4), xvar(-3), xvar(-2), xvar(-1), xvar(0), xvar(1), xvar(2), xvar(3), xvar(4));
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = interp_10(yvar(-5), yvar(-4), yvar(-3), yvar(-2), yvar(-1), yvar(0), yvar(1), yvar(2), yvar(3), yvar(4));
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = interp_10(zvar(-5), zvar(-4), zvar(-3), zvar(-2), zvar(-1), zvar(0), zvar(1), zvar(2), zvar(3), zvar(4));
      }
    }
  });
}

#endif
