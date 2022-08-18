
#pragma once

#include "cfv_recon.h"
#include "common.h"
#include "weno_func_recon.h"
#include "weno_recon.h"

template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord = 2,
          uint hs = (ord - 1) / 2>
void YAKL_INLINE compute_twisted_edge_recon(
    const real5d &edgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 3, ord, ord, ord> const &wenoRecon,
    SArray<real, 2, ord, tord> const &to_gll,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {

  SArray<real, 3, ndofs, ndims, ord> stencil;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          stencil(l, d, p) = var(l, k + ks, j + js, i + is + p - hs, n);
        }
        if (d == 1) {
          stencil(l, d, p) = var(l, k + ks, j + js + p - hs, i + is, n);
        }
      }
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv<ndofs, ndims>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    weno<ndofs, ndims>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func<ndofs, ndims, ord>(edgerecon, stencil, wenoRecon, to_gll, wenoIdl,
                                 wenoSigma);
  }

  for (int d = 0; d < ndims; d++) {
    for (int l = 0; l < ndofs; l++) {
      for (int m = 0; m < 2; m++) {
        edgereconvar(l + d * ndofs + ndofs * ndims * m, k + ks, j + js, i + is,
                     n) = edgerecon(l, d, m);
      }
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord = 2,
          uint hs = (ord - 1) / 2>
void YAKL_INLINE compute_twisted_vert_edge_recon(
    const real5d &vertedgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 3, ord, ord, ord> const &wenoRecon,
    SArray<real, 2, ord, tord> const &to_gll,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {
  SArray<real, 3, ndofs, 1, ord> stencil;
  SArray<real, 3, ndofs, 1, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      stencil(l, 0, p) = var(l, k + ks + p - hs, j + js, i + is, n);
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv<ndofs, 1>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    weno<ndofs, 1>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func<ndofs, 1, ord>(edgerecon, stencil, wenoRecon, to_gll, wenoIdl,
                             wenoSigma);
  }

  for (int l = 0; l < ndofs; l++) {
    for (int m = 0; m < 2; m++) {
      vertedgereconvar(l + ndofs * m, k + ks, j + js, i + is, n) =
          edgerecon(l, 0, m);
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord = 2,
          uint hs = (ord - 1) / 2>
void YAKL_INLINE compute_straight_edge_recon(
    const real5d &edgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 3, ord, ord, ord> const &wenoRecon,
    SArray<real, 2, ord, tord> const &to_gll,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {

  SArray<real, 3, ndofs, ndims, ord> stencil;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = ndims - 1; d >= 0; d--) {
        if (d == ndims - 1) {
          stencil(l, d, p) = var(l, k + ks, j + js, i + is + p - hs, n);
        }
        if (d == ndims - 2) {
          stencil(l, d, p) = var(l, k + ks, j + js + p - hs, i + is, n);
        }
        // if (d==ndims-3) {stencil(l,d,p) = var(l, k+ks+p-hs, j+js, i+is);}
      }
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv<ndofs, ndims>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    weno<ndofs, ndims>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func<ndofs, ndims, ord>(edgerecon, stencil, wenoRecon, to_gll, wenoIdl,
                                 wenoSigma);
  }

  for (int d = ndims - 1; d >= 0; d--) {
    for (int l = 0; l < ndofs; l++) {
      for (int m = 0; m < 2; m++) {
        edgereconvar(l + d * ndofs + ndofs * ndims * m, k + ks, j + js, i + is,
                     n) = edgerecon(l, d, m);
      }
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord = 2,
          uint hs = (ord - 1) / 2>
void YAKL_INLINE compute_straight_xz_edge_recon(
    const real5d &edgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 3, ord, ord, ord> const &wenoRecon,
    SArray<real, 2, ord, tord> const &to_gll,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {
  SArray<real, 3, ndofs, 1, ord> stencil;
  SArray<real, 3, ndofs, 1, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      // The +1 in k here is required since twisted 0-forms have extra dofs at
      // the top and bottom, and we are looping over straight cells!
      stencil(l, 0, p) = var(l, k + ks + 1, j + js, i + is + p - hs, n);
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv<ndofs, 1>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    weno<ndofs, 1>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func<ndofs, 1, ord>(edgerecon, stencil, wenoRecon, to_gll, wenoIdl,
                             wenoSigma);
  }

  for (int l = 0; l < ndofs; l++) {
    for (int m = 0; m < 2; m++) {
      edgereconvar(l + ndofs * m, k + ks, j + js, i + is, n) =
          edgerecon(l, 0, m);
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype, uint ord, uint tord = 2,
          uint hs = (ord - 1) / 2>
void YAKL_INLINE compute_straight_xz_vert_edge_recon(
    const real5d &edgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 3, ord, ord, ord> const &wenoRecon,
    SArray<real, 2, ord, tord> const &to_gll,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {
  SArray<real, 3, ndofs, 1, ord> stencil;
  SArray<real, 3, ndofs, 1, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      // The +1 in k here is required since twisted 0-forms have extra dofs at
      // the top and bottom, and we are looping over straight cells!
      stencil(l, 0, p) = var(l, k + ks + p - hs + 1, j + js, i + is, n);
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv<ndofs, 1>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    weno<ndofs, 1>(edgerecon, stencil);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func<ndofs, 1, ord>(edgerecon, stencil, wenoRecon, to_gll, wenoIdl,
                             wenoSigma);
  }

  for (int l = 0; l < ndofs; l++) {
    for (int m = 0; m < 2; m++) {
      edgereconvar(l + ndofs * m, k + ks, j + js, i + is, n) =
          edgerecon(l, 0, m);
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE
centered_recon(SArray<real, 2, ndofs, nd> &recon,
               SArray<real, 3, ndofs, nd, 2> const &edgerecon) {

  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      recon(l, d) = 0.5 * (edgerecon(l, d, 1) + edgerecon(l, d, 0));
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE upwind_recon(SArray<real, 2, ndofs, nd> &recon,
                              SArray<real, 3, ndofs, nd, 2> const &edgerecon,
                              SArray<real, 1, nd> const &flux) {

  real upwind_param;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      upwind_param = copysign(1.0, flux(d));
      upwind_param = 0.5 * (upwind_param + fabs(upwind_param));
      recon(l, d) = edgerecon(l, d, 1) * (1. - upwind_param) +
                    edgerecon(l, d, 0) * upwind_param;
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_twisted_recon(const real5d &reconvar,
                                       const real5d &edgereconvar,
                                       const real5d &U, int is, int js, int ks,
                                       int i, int j, int k, int n) {

  SArray<real, 2, ndofs, ndims> recon;
  SArray<real, 1, ndims> uvar;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int d = 0; d < ndims; d++) {
    uvar(d) = U(d, k + ks, j + js, i + is, n);
    for (int l = 0; l < ndofs; l++) {
      if (d == 0) {
        edgerecon(l, d, 0) = edgereconvar(l + d * ndofs + ndofs * ndims * 1,
                                          k + ks, j + js, i + is - 1, n);
        edgerecon(l, d, 1) = edgereconvar(l + d * ndofs + ndofs * ndims * 0,
                                          k + ks, j + js, i + is, n);
      }
      if (d == 1) {
        edgerecon(l, d, 0) = edgereconvar(l + d * ndofs + ndofs * ndims * 1,
                                          k + ks, j + js - 1, i + is, n);
        edgerecon(l, d, 1) = edgereconvar(l + d * ndofs + ndofs * ndims * 0,
                                          k + ks, j + js, i + is, n);
      }
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, ndims>(recon, edgerecon);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
  }

  for (int d = 0; d < ndims; d++) {
    for (int l = 0; l < ndofs; l++) {
      reconvar(l + d * ndofs, k + ks, j + js, i + is, n) = recon(l, d);
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_twisted_vert_recon(const real5d &vertreconvar,
                                            const real5d &vertedgereconvar,
                                            const real5d &UW, int is, int js,
                                            int ks, int i, int j, int k,
                                            int n) {

  SArray<real, 2, ndofs, 1> recon;
  SArray<real, 1, 1> uvar;
  SArray<real, 3, ndofs, 1, 2> edgerecon;

  uvar(0) = UW(0, k + ks, j + js, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    edgerecon(l, 0, 0) =
        vertedgereconvar(l + ndofs * 1, k + ks - 1, j + js, i + is, n);
    edgerecon(l, 0, 1) =
        vertedgereconvar(l + ndofs * 0, k + ks, j + js, i + is, n);
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, 1>(recon, edgerecon);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
  }

  for (int l = 0; l < ndofs; l++) {
    vertreconvar(l, k + ks, j + js, i + is, n) = recon(l, 0);
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_straight_recon(const real5d &reconvar,
                                        const real5d &edgereconvar,
                                        const real5d &UT, int is, int js,
                                        int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndofs, ndims> recon;
  SArray<real, 1, ndims> uvar;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int d = ndims - 1; d >= 0; d--) {
    uvar(d) = UT(d, k + ks, j + js, i + is, n);
    if (ndims == 2 && d == ndims - 2) {
      uvar(d) = -uvar(d);
    } // corrects for "twist" in 2D

    for (int l = 0; l < ndofs; l++) {
      if (d == ndims - 1) {
        edgerecon(l, d, 0) = edgereconvar(l + d * ndofs + ndofs * ndims * 1,
                                          k + ks, j + js, i + is, n);
        edgerecon(l, d, 1) = edgereconvar(l + d * ndofs + ndofs * ndims * 0,
                                          k + ks, j + js, i + is + 1, n);
      }
      if (d == ndims - 2) {
        edgerecon(l, d, 0) = edgereconvar(l + d * ndofs + ndofs * ndims * 1,
                                          k + ks, j + js, i + is, n);
        edgerecon(l, d, 1) = edgereconvar(l + d * ndofs + ndofs * ndims * 0,
                                          k + ks, j + js + 1, i + is, n);
      }
      // if (d==ndims-3)
      // {
      // edgerecon(l,d,0) = edgereconvar(l + d*ndofs + ndofs*ndims*1, k+ks,
      // j+js, i+is); edgerecon(l,d,1) = edgereconvar(l + d*ndofs +
      // ndofs*ndims*0, k+ks+1, j+js, i+is);
      // }
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, ndims>(recon, edgerecon);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
  }

  for (int d = ndims - 1; d >= 0; d--) {
    for (int l = 0; l < ndofs; l++) {
      reconvar(l + d * ndofs, k + ks, j + js, i + is, n) = recon(l, d);
    }
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_straight_xz_recon(const real5d &reconvar,
                                           const real5d &edgereconvar,
                                           const real5d &WT, int is, int js,
                                           int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndofs, 1> recon;
  SArray<real, 1, 1> uvar;
  SArray<real, 3, ndofs, 1, 2> edgerecon;
  uvar(0) = WT(0, k + ks, j + js, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    edgerecon(l, 0, 0) = edgereconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    edgerecon(l, 0, 1) =
        edgereconvar(l + ndofs * 0, k + ks, j + js, i + is + 1, n);
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, 1>(recon, edgerecon);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
  }

  for (int l = 0; l < ndofs; l++) {
    reconvar(l, k + ks, j + js, i + is, n) = recon(l, 0);
  }
}

template <uint ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_straight_xz_vert_recon(const real5d &reconvar,
                                                const real5d &edgereconvar,
                                                const real5d &VT, int is,
                                                int js, int ks, int i, int j,
                                                int k, int n) {
  SArray<real, 2, ndofs, 1> recon;
  SArray<real, 1, 1> uvar;
  SArray<real, 3, ndofs, 1, 2> edgerecon;
  // Needs a "twist"
  uvar(0) = -VT(0, k + ks, j + js, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    edgerecon(l, 0, 0) =
        edgereconvar(l + ndofs * 1, k + ks - 1, j + js, i + is, n);
    edgerecon(l, 0, 1) = edgereconvar(l + ndofs * 0, k + ks, j + js, i + is, n);
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, 1>(recon, edgerecon);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENO) {
    upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
  }
  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
  }

  for (int l = 0; l < ndofs; l++) {
    reconvar(l, k + ks, j + js, i + is, n) = recon(l, 0);
  }
}
