
#pragma once

#include "cfv_recon.h"
#include "common.h"
#include "weno_func_recon.h"
#include "weno_func_recon_variable.h"
#include "weno_recon.h"

namespace pamc {

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t tord = 2, index_t hs = (ord - 1) / 2>
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

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t tord = 2, index_t hs = (ord - 1) / 2>
void YAKL_INLINE compute_twisted_vert_edge_recon_uniform(
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

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t hs = (ord - 1) / 2>
void YAKL_INLINE compute_twisted_vert_edge_recon_variable(
    const real5d &vertedgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 2, ord, 2> const &coefs_to_gll,
    SArray<real, 2, ord, 2> const &sten_to_gll,
    SArray<real, 2, ord, ord> const &sten_to_coefs,
    SArray<real, 3, hs + 1, hs + 1, hs + 1> const &weno_recon_lower,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {
  SArray<real, 2, ndofs, ord> stencil;
  SArray<real, 2, ndofs, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      stencil(l, p) = var(l, k + ks + p - hs, j + js, i + is, n);
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv_vert<ndofs, ord>(edgerecon, stencil);
  }

  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func_vert<ndofs, ord>(edgerecon, stencil, coefs_to_gll, sten_to_gll,
                               sten_to_coefs, weno_recon_lower, wenoIdl,
                               wenoSigma);
  }

  for (int l = 0; l < ndofs; l++) {
    for (int m = 0; m < 2; m++) {
      vertedgereconvar(l + ndofs * m, k + ks, j + js, i + is, n) =
          edgerecon(l, m);
    }
  }
}

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t tord = 2, index_t hs = (ord - 1) / 2>
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

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t tord = 2, index_t hs = (ord - 1) / 2>
void YAKL_INLINE compute_straight_hz_edge_recon(
    const real5d &edgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 3, ord, ord, ord> const &wenoRecon,
    SArray<real, 2, ord, tord> const &to_gll,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {
  SArray<real, 3, ndofs, ndims, ord> stencil;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        // The +1 in k here is required since twisted 0-forms have extra dofs at
        // the top and bottom, and we are looping over straight cells!
        if (d == 0) {
          stencil(l, d, p) =
              var(l + ndofs * d, k + ks + 1, j + js, i + is + p - hs, n);
        }
        if (d == 1) {
          stencil(l, d, p) =
              var(l + ndofs * d, k + ks + 1, j + js + p - hs, i + is, n);
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
        edgereconvar(l + ndofs * d + ndofs * ndims * m, k + ks, j + js, i + is,
                     n) = edgerecon(l, d, m);
      }
    }
  }
}

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t tord = 2, index_t hs = (ord - 1) / 2>
void YAKL_INLINE compute_straight_hz_vert_edge_recon_uniform(
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

template <index_t ndofs, RECONSTRUCTION_TYPE recontype, index_t ord,
          index_t hs = (ord - 1) / 2>
void YAKL_INLINE compute_straight_hz_vert_edge_recon_variable(
    const real5d &edgereconvar, const real5d &var, int is, int js, int ks,
    int i, int j, int k, int n, SArray<real, 2, ord, 2> const &coefs_to_gll,
    SArray<real, 2, ord, 2> const &sten_to_gll,
    SArray<real, 2, ord, ord> const &sten_to_coefs,
    SArray<real, 3, hs + 1, hs + 1, hs + 1> const &weno_recon_lower,
    SArray<real, 1, hs + 2> const &wenoIdl, real wenoSigma) {
  SArray<real, 2, ndofs, ord> stencil;
  SArray<real, 2, ndofs, 2> edgerecon;

  for (int p = 0; p < ord; p++) {
    for (int l = 0; l < ndofs; l++) {
      // The +1 in k here is required since twisted 0-forms have extra dofs at
      // the top and bottom, and we are looping over straight cells!
      stencil(l, p) = var(l, k + ks + p - hs + 1, j + js, i + is, n);
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    cfv_vert<ndofs, ord>(edgerecon, stencil);
  }

  if (recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    weno_func_vert<ndofs, ord>(edgerecon, stencil, coefs_to_gll, sten_to_gll,
                               sten_to_coefs, weno_recon_lower, wenoIdl,
                               wenoSigma);
  }

  for (int l = 0; l < ndofs; l++) {
    for (int m = 0; m < 2; m++) {
      edgereconvar(l + ndofs * m, k + ks, j + js, i + is, n) = edgerecon(l, m);
    }
  }
}

template <index_t ndofs, index_t nd>
void YAKL_INLINE
centered_recon(SArray<real, 2, ndofs, nd> &recon,
               SArray<real, 3, ndofs, nd, 2> const &edgerecon) {

  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      recon(l, d) = 0.5 * (edgerecon(l, d, 1) + edgerecon(l, d, 0));
    }
  }
}

template <index_t ndofs, index_t nd>
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

template <index_t ndofs, index_t nd>
void YAKL_INLINE
tanh_upwind_recon(SArray<real, 2, ndofs, nd> &recon,
                  SArray<real, 3, ndofs, nd, 2> const &edgerecon,
                  SArray<real, 1, nd> const &flux, real tanh_upwind_coeff) {

  real upwind_param;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      upwind_param = std::tanh(flux(d) * tanh_upwind_coeff);
      recon(l, d) = 0.5_fp * (edgerecon(l, d, 1) * (1 - upwind_param) +
                              edgerecon(l, d, 0) * (upwind_param + 1));
    }
  }
}

template <index_t ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_twisted_recon(const real5d &reconvar,
                                       const real5d &edgereconvar,
                                       const Geometry<Twisted> &dgeom,
                                       const real5d &Fvar,
                                       real tanh_upwind_coeff, int is, int js,
                                       int ks, int i, int j, int k, int n) {

  SArray<real, 2, ndofs, ndims> recon;
  SArray<real, 1, ndims> uvar;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int d = 0; d < ndims; d++) {
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
  } else if (recontype == RECONSTRUCTION_TYPE::WENO ||
             recontype == RECONSTRUCTION_TYPE::WENOFUNC) {

    for (int d = 0; d < ndims; ++d) {
      uvar(d) = Fvar(d, k + ks, j + js, i + is, n);
    }

    if (dual_upwind_type == UPWIND_TYPE::HEAVISIDE) {
      upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
    } else if (dual_upwind_type == UPWIND_TYPE::TANH) {
      for (int d = 0; d < ndims; ++d) {
        uvar(d) /= dgeom.get_area_nm11entity(d, k + ks, j + js, i + is, n);
      }
      tanh_upwind_recon<ndofs, ndims>(recon, edgerecon, uvar,
                                      tanh_upwind_coeff);
    }
  }

  for (int d = 0; d < ndims; d++) {
    for (int l = 0; l < ndofs; l++) {
      reconvar(d + l * ndims, k + ks, j + js, i + is, n) = recon(l, d);
    }
  }
}

template <index_t ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_twisted_vert_recon(
    const real5d &vertreconvar, const real5d &vertedgereconvar,
    const Geometry<Twisted> &dgeom, const real5d &FWvar, real tanh_upwind_coeff,
    int is, int js, int ks, int i, int j, int k, int n) {

  SArray<real, 2, ndofs, 1> recon;
  SArray<real, 1, 1> uvar;
  SArray<real, 3, ndofs, 1, 2> edgerecon;

  for (int l = 0; l < ndofs; l++) {
    edgerecon(l, 0, 0) =
        vertedgereconvar(l + ndofs * 1, k + ks - 1, j + js, i + is, n);
    edgerecon(l, 0, 1) =
        vertedgereconvar(l + ndofs * 0, k + ks, j + js, i + is, n);
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, 1>(recon, edgerecon);
  } else if (recontype == RECONSTRUCTION_TYPE::WENO ||
             recontype == RECONSTRUCTION_TYPE::WENOFUNC) {

    uvar(0) = FWvar(0, k + ks, j + js, i + is, n);

    if (dual_vert_upwind_type == UPWIND_TYPE::HEAVISIDE) {
      upwind_recon<ndofs, 1>(recon, edgerecon, uvar);
    } else if (dual_vert_upwind_type == UPWIND_TYPE::TANH) {
      uvar(0) /= dgeom.get_area_n0entity(k + ks, j + js, i + is, n);
      tanh_upwind_recon<ndofs, 1>(recon, edgerecon, uvar, tanh_upwind_coeff);
    }
  }

  for (int l = 0; l < ndofs; l++) {
    vertreconvar(l, k + ks, j + js, i + is, n) = recon(l, 0);
  }
}

template <index_t ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_straight_recon(const real5d &reconvar,
                                        const real5d &edgereconvar,
                                        const Geometry<Straight> &pgeom,
                                        const real5d &UT,
                                        real tanh_upwind_coeff, int is, int js,
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
  } else if (recontype == RECONSTRUCTION_TYPE::WENO ||
             recontype == RECONSTRUCTION_TYPE::WENOFUNC) {
    if (upwind_type == UPWIND_TYPE::HEAVISIDE) {
      upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
    } else {
      for (int d = 0; d < ndims; ++d) {
        uvar(d) /= pgeom.get_area_nm11entity(d, k + ks, j + js, is, n);
      }
      tanh_upwind_recon<ndofs, ndims>(recon, edgerecon, uvar,
                                      tanh_upwind_coeff);
    }
  }

  for (int d = ndims - 1; d >= 0; d--) {
    for (int l = 0; l < ndofs; l++) {
      reconvar(l + d * ndofs, k + ks, j + js, i + is, n) = recon(l, d);
    }
  }
}

#ifdef PAMC_EXTRUDED
template <index_t ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_straight_hz_recon(const real5d &reconvar,
                                           const real5d &edgereconvar,
                                           const Geometry<Straight> &pgeom,
                                           const real5d &FTWvar,
                                           real tanh_upwind_coeff, int is,
                                           int js, int ks, int i, int j, int k,
                                           int n) {
  SArray<real, 2, ndofs, ndims> recon;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int d = 0; d < ndims; d++) {
    for (int l = 0; l < ndofs; l++) {
      if (d == 0) {
        edgerecon(l, d, 0) = edgereconvar(l + ndofs * d + ndofs * ndims * 1,
                                          k + ks, j + js, i + is, n);
        edgerecon(l, d, 1) = edgereconvar(l + ndofs * d + ndofs * ndims * 0,
                                          k + ks, j + js, i + is + 1, n);
      }
      if (d == 1) {
        edgerecon(l, d, 0) = edgereconvar(l + ndofs * d + ndofs * ndims * 1,
                                          k + ks, j + js, i + is, n);
        edgerecon(l, d, 1) = edgereconvar(l + ndofs * d + ndofs * ndims * 0,
                                          k + ks, j + js + 1, i + is, n);
      }
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, ndims>(recon, edgerecon);
  } else if (recontype == RECONSTRUCTION_TYPE::WENO ||
             recontype == RECONSTRUCTION_TYPE::WENOFUNC) {

    SArray<real, 1, ndims> uvar;
    for (int d = 0; d < ndims; ++d) {
      uvar(d) = FTWvar(d, k + ks, j + js, i + is, n);
    }

    if (upwind_type == UPWIND_TYPE::HEAVISIDE) {
      upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
    } else if (upwind_type == UPWIND_TYPE::TANH) {
      for (int d = 0; d < ndims; ++d) {
        uvar(d) /= pgeom.get_area_nm11entity(d, k + ks, j + js, is, n);
      }
      tanh_upwind_recon<ndofs, ndims>(recon, edgerecon, uvar,
                                      tanh_upwind_coeff);
    }
  }

  for (int d = 0; d < ndims; d++) {
    for (int l = 0; l < ndofs; l++) {
      reconvar(l + ndofs * d, k + ks, j + js, i + is, n) = recon(l, d);
    }
  }
}

template <index_t ndofs, RECONSTRUCTION_TYPE recontype>
void YAKL_INLINE compute_straight_hz_vert_recon(const real5d &reconvar,
                                                const real5d &edgereconvar,
                                                const Geometry<Straight> &pgeom,
                                                const real5d &FTvar,
                                                real tanh_upwind_coeff, int is,
                                                int js, int ks, int i, int j,
                                                int k, int n) {
  SArray<real, 2, ndofs, ndims> recon;
  SArray<real, 3, ndofs, ndims, 2> edgerecon;

  for (int d = 0; d < ndims; d++) {
    for (int l = 0; l < ndofs; l++) {
      edgerecon(l, d, 0) = edgereconvar(l + ndofs * d + ndofs * ndims * 1,
                                        k + ks - 1, j + js, i + is, n);
      edgerecon(l, d, 1) = edgereconvar(l + ndofs * d + ndofs * ndims * 0,
                                        k + ks, j + js, i + is, n);
    }
  }

  if (recontype == RECONSTRUCTION_TYPE::CFV) {
    centered_recon<ndofs, ndims>(recon, edgerecon);
  } else if (recontype == RECONSTRUCTION_TYPE::WENO ||
             recontype == RECONSTRUCTION_TYPE::WENOFUNC) {

    SArray<real, 1, ndims> uvar;
    for (int d = 0; d < ndims; ++d) {
      uvar(d) = FTvar(d, k + ks, j + js, i + is, n);
      if (d == 0) {
        uvar(d) *= -1; // needs a "twist"
      }
    }

    if (vert_upwind_type == UPWIND_TYPE::HEAVISIDE) {
      upwind_recon<ndofs, ndims>(recon, edgerecon, uvar);
    } else if (vert_upwind_type == UPWIND_TYPE::TANH) {
      uvar(0) /= pgeom.get_area_n0entity(k + ks, j + js, i + is, n);
      tanh_upwind_recon<ndofs, ndims>(recon, edgerecon, uvar,
                                      tanh_upwind_coeff);
    }
  }

  for (int d = 0; d < ndims; ++d) {
    for (int l = 0; l < ndofs; l++) {
      reconvar(l + ndofs * d, k + ks, j + js, i + is, n) = recon(l, d);
    }
  }
}
#endif
} // namespace pamc
