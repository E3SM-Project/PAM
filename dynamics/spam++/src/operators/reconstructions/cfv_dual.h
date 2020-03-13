
#ifndef _CFV_H_
#define _CFV_H_

#include "common.h"
#include "topology.h"
#include "geometry.h"

//ALL OF THESE ASSUME A DIAGONAL HODGE STAR THAT MAPS FROM N-FORMS TO 0-FORMS
//WE LIKELY WANT TO BE MORE CLEVER FOR HIGHER ORDER VERSIONS...



#define xvar(p) var(l, k+ks, j+js, i+is+p)/geom.get_J_cell(k+ks, j+js, i+is+p)
#define yvar(p) var(l, k+ks, j+js+p, i+is)/geom.get_J_cell(k+ks, j+js+p, i+is)
#define zvar(p) var(l, k+ks+p, j+js, i+is)/geom.get_J_cell(k+ks+p, j+js, i+is)




// CFV Interpolations
real YAKL_INLINE interp_2(real phi, real phip1){
    return 0.5*(phi + phip1);
};

real YAKL_INLINE interp_4(real phim1, real phi, real phip1, real phip2){
    return (7.0/12.0)*(phi + phip1 ) -(1.0/12.0)*(phim1 + phip2);
};


real YAKL_INLINE interp_6(real phim2, real phim1, real phi, real phip1,
                real phip2, real phip3){
    return ((37.0/60.0) * (phi + phip1) - (2.0/15.0)*(phim1 + phip2)
            + (1.0/60.0)*(phim2 + phip3));
};


real YAKL_INLINE interp_8(real phim3, real phim2, real phim1, real phi,
                real phip1, real phip2, real phip3, real phip4){
   return  (533./840. * (phi + phip1) - 139.0/840.0 * (phim1 + phip2 )
            + 29.0/840.0 * (phim2 + phip3) -1.0/280.0*(phim3 + phip4));
};


real YAKL_INLINE interp_10(real phim4, real phim3, real phim2, real phim1,
                real phi, real phip1, real phip2, real phip3,
                real phip4, real phip5){
    return (1627.0/2520.0 * (phi + phip1) - 473.0/2520.0 * (phim1 + phip2 )
                       + 127.0/2520.0* (phim2 + phip3) -23.0/2520.0 *(phim3 + phip4)
                        + 1.0/1260.0*(phim4 + phip5));
};













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
