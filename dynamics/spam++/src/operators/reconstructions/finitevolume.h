
#ifndef _FINITEVOLUME_H_
#define _FINITEVOLUME_H_

#include "common.h"
#include "topology.h"
#include "geometry.h"

//THIS REALLY ASSUMES A DIAGONAL HODGE STAR THAT MAPS FROM N-FORMS TO 0-FORMS
//WE LIKELY WANT TO BE MORE CLEVER FOR HIGHER ORDER VERSIONS...

template<uint ndims, uint ndofs> void YAKL_INLINE ufv1_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  real upwind_param;
  yakl::parallel_for("ComputeUFV1Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      upwind_param = copysign(1.0, flux(0, k+ks, j+js, i+is));
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) * (1. - upwind_param) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1) * upwind_param);

      //y-dir
      if (ndims >= 2) {
      upwind_param = copysign(1.0, flux(1, k+ks, j+js, i+is));
      recon(l+1*ndofs, k+ks, j+js, i+is) =  (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) * (1. - upwind_param) + var(l, k+ks,   j+js-1,   i+is)/geom.get_J_cell(k+ks, j+js-1, i+is) * upwind_param);
      }
      //z-dir
      if (ndims == 3) {
      upwind_param = copysign(1.0, flux(2, k+ks, j+js, i+is));
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) * (1. - upwind_param) + var(l, k+ks-1,   j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is) * upwind_param);
      }
    }
  });
}





















template<uint ndims, uint ndofs> void YAKL_INLINE cfv2_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV2Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))/2.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))/2.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))/2.0;
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
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}

//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE cfv6_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV6Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;

      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}

//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE cfv8_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV8Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}


//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE cfv10_recon(realArr recon, const realArr var, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeCFV10Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}



























//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE weno3_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeWENO3Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}


//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE weno5_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeWENO5Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}


//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE weno7_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeWENO7Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}

//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE weno9_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeWENO9Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}


//FIX THIS
template<uint ndims, uint ndofs> void YAKL_INLINE weno11_recon(realArr recon, const realArr var, const realArr flux, const Topology<ndims> &topology, Geometry<ndims,1,1,1> &geom) {

  int is = topology.is;
  int js = topology.js;
  int ks = topology.ks;

  yakl::parallel_for("ComputeWENO11Recon", topology.n_cells, YAKL_LAMBDA (int iGlob) {
    int k, j, i;
    yakl::unpackIndices(iGlob, topology.n_cells_z, topology.n_cells_y, topology.n_cells_x, k, j, i);
    for (int l=0; l<ndofs; l++) {

      //x-dir
      recon(l+0*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js,   i+is-1)/geom.get_J_cell(k+ks, j+js, i+is-1))*7./12.0
                                          +(var(l, k+ks, j+js, i+is+1)/geom.get_J_cell(k+ks, j+js, i+is+1) + var(l, k+ks,   j+js,   i+is-2)/geom.get_J_cell(k+ks, j+js, i+is-2))*1./12.0;
      //y-dir
      if (ndims >= 2) {
      recon(l+1*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks,   j+js-1, i+is)/geom.get_J_cell(k+ks, j+js-1, i+is))*7./12.0;
                                          +(var(l, k+ks, j+js+1, i+is)/geom.get_J_cell(k+ks, j+js+1, i+is) + var(l, k+ks,   j+js-2,   i+is)/geom.get_J_cell(k+ks, j+js-2, i+is))*1./12.0;
      }
      //z-dir
      if (ndims == 3) {
      recon(l+2*ndofs, k+ks, j+js, i+is) = (var(l, k+ks, j+js, i+is)/geom.get_J_cell(k+ks, j+js, i+is) + var(l, k+ks-1, j+js,   i+is)/geom.get_J_cell(k+ks-1, j+js, i+is))*7./12.0;
                                          +(var(l, k+ks+1, j+js, i+is)/geom.get_J_cell(k+ks+1, j+js, i+is) + var(l, k+ks-2,   j+js,   i+is)/geom.get_J_cell(k+ks-2, j+js, i+is))*1./12.0;
      }
    }
  });
}


#endif
