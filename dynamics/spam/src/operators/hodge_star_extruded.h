#pragma once

#include "common.h"
#include "geometry.h"
#include "hodge_star.h"

namespace pamc {

real YAKL_INLINE H00_diagonal(const Geometry<Straight> &pgeom,
                              const Geometry<Twisted> &dgeom, int is, int js,
                              int ks, int i, int j, int k, int n) {
  return dgeom.get_area_n1entity(k + ks, j + js, i + is, n) /
         pgeom.get_area_00entity(k + ks, j + js, i + is, n);
}

template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H00(SArray<real, 1, ndofs> &x0, const real5d &var,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real H00geom = H00_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int p = 0; p < hord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks, j + js, i + is + p - hoff, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks, j + js + p - hoff, i + is, n);
        }
      }
    }
  }
  H0<ndofs>(x0, x, H00geom);
}

template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H00(const real5d &var0, const real5d &var,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H00<ndofs, hord, vord, hoff, voff>(x0, var, pgeom, dgeom, is, js, ks,
                                             i, j, k, n);

  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}

real YAKL_INLINE H00bar_diagonal(const Geometry<Straight> &pgeom,
                                 const Geometry<Twisted> &dgeom, int is, int js,
                                 int ks, int i, int j, int k, int n) {
  return pgeom.get_area_n1entity(k + ks, j + js, i + is, n) /
         dgeom.get_area_00entity(k + ks, j + js, i + is, n);
}

template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H00bar(SArray<real, 1, ndofs> &x0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real H00bargeom = H00bar_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int p = 0; p < hord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks, j + js, i + is + p - hoff, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks, j + js + p - hoff, i + is, n);
        }
      }
    }
  }
  H0bar<ndofs>(x0, x, H00bargeom);
}

template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H00bar(const real5d &var0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H00bar<ndofs, hord, vord, hoff, voff>(x0, var, pgeom, dgeom, is, js,
                                                ks, i, j, k, n);

  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}

void YAKL_INLINE H10_diagonal(SArray<real, 1, ndims> &H10diag,
                              const Geometry<Straight> &pgeom,
                              const Geometry<Twisted> &dgeom, int is, int js,
                              int ks, int i, int j, int k, int n) {
  for (int d = 0; d < ndims; d++) {
    H10diag(d) = dgeom.get_area_nm11entity(d, k + ks, j + js, i + is, n) /
                 pgeom.get_area_10entity(d, k + ks, j + js, i + is, n);
  }
}

// BROKEN FOR 2D+1D EXT
// MAINLY IN THE AREA CALCS...
template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H10(SArray<real, 2, ndofs, ndims> &u,
                             const real5d &vvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> v;

  SArray<real, 1, ndims> H10geom;
  H10_diagonal(H10geom, pgeom, dgeom, is, js, ks, i, j, k, n);

  for (int l = 0; l < ndofs; l++) {
    for (int p = 0; p < ord - 1; p++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          v(l, d, p) = vvar(l * ndims + d, k + ks, j + js, i + is + p - off, n);
        }
        if (d == 1) {
          v(l, d, p) = vvar(l * ndims + d, k + ks, j + js + p - off, i + is, n);
        }
      }
    }
  }
  H1(u, v, H10geom);
}

// No indexing issues since we go from p10 to d01, and d01 has no "extended
// boundary"
template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H10(const real5d &uvar, const real5d &vvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndofs, ndims> u;
  compute_H10<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {
        uvar(l * ndims + d, k + ks, j + js, i + is, n) = u(l, d);
      }
      if (addmode == ADD_MODE::ADD) {
        uvar(l * ndims + d, k + ks, j + js, i + is, n) += u(l, d);
      }
    }
  }
}

template <uint ord, int off = ord / 2 - 1>
void YAKL_INLINE fourier_H10(SArray<real, 1, ndims> &u,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n, int nx, int ny,
                             int nz) {

  SArray<real, 2, ndims, off + GPU_PAD> shift;

  SArray<real, 1, ndims> H10geom;
  H10_diagonal(H10geom, pgeom, dgeom, is, js, ks, i, j, k, n);

  for (int p = 0; p < off; p++) {
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        shift(d, p) = (2 * pi * i * (p - off)) / nx;
      }
      if (d == 1) {
        shift(d, p) = (2 * pi * j * (p - off)) / ny;
      }
    }
  }
  H1hat(u, H10geom, shift);
}

// Note the indexing here, this is key
real YAKL_INLINE H01_diagonal(const Geometry<Straight> &pgeom,
                              const Geometry<Twisted> &dgeom, int is, int js,
                              int ks, int i, int j, int k, int n) {
  // THIS IS 2ND ORDER AT BEST...
  return dgeom.get_area_n0entity(k + ks, j + js, i + is, n) /
         pgeom.get_area_01entity(k + ks - 1, j + js, i + is, n);
}

template <uint ndofs>
void YAKL_INLINE compute_H01(SArray<real, 1, ndofs> &uw, const real5d &wvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  for (int l = 0; l < ndofs; ++l) {
    const real H01_geom = H01_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
    uw(l) = wvar(l, k + ks - 1, j + js, i + is, n) * H01_geom;
  }
}

// Indexing issues since we go from p01 to d10, and d10 has an "extended
// boundary" ie boundary vert edges Since we index over d10, need to subtract 1
// from k when indexing p01 ie the kth edge flux corresponds with the k-1th edge
// velocity Also should be called with k=[1,...,ni-2] ie skip the first and last
// fluxes, which are set diagnostically (=0 for no-flux bcs)
template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H01(const real5d &uwvar, const real5d &wvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> uw;
  compute_H01(uw, wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      uwvar(l, k + ks, j + js, i + is, n) = uw(l);
    }
    if (addmode == ADD_MODE::ADD) {
      uwvar(l, k + ks, j + js, i + is, n) += uw(l);
    }
  }
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H01(SArray<real, 1, 1> &uwvar, const real5d &wvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> uw;
  compute_H01(uw, wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      uwvar(l) = uw(l);
    }
    if (addmode == ADD_MODE::ADD) {
      uwvar(l) += uw(l);
    }
  }
}

void YAKL_INLINE Hnm11_diagonal(SArray<real, 1, ndims> &Hnm11diag,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  for (int d = 0; d < ndims; d++) {
    Hnm11diag(d) = dgeom.get_area_10entity(d, k + ks, j + js, i + is, n) /
                   pgeom.get_area_nm11entity(d, k - 1 + ks, j + js, i + is, n);
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hnm11(SArray<real, 1, ndims> &u, const real5d &vvar,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;

  SArray<real, 1, ndims> Hnm11geom;
  Hnm11_diagonal(Hnm11geom, pgeom, dgeom, is, js, ks, i, j, k, n);

  for (int p = 0; p < ord - 1; p++) {
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        v(d, p) = vvar(d, k - 1 + ks, j + js, i + is + p - off, n);
      }
      if (d == 1) {
        v(d, p) = vvar(d, k - 1 + ks, j + js + p - off, i + is, n);
      }
    }
  }
  H1bar(u, v, Hnm11geom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hnm11(const real5d &uvar, const real5d &vvar,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> u;
  compute_Hnm11<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      uvar(d, k + ks, j + js, i + is, n) = u(d);
    }
    if (addmode == ADD_MODE::ADD) {
      uvar(d, k + ks, j + js, i + is, n) += u(d);
    }
  }
}

real YAKL_INLINE Hn0_diagonal(const Geometry<Straight> &pgeom,
                              const Geometry<Twisted> &dgeom, int is, int js,
                              int ks, int i, int j, int k, int n) {
  return -dgeom.get_area_01entity(k + ks, j + js, i + is, n) /
         pgeom.get_area_n0entity(k + ks, j + js, i + is, n);
}

real YAKL_INLINE compute_Hn0(const real5d &wvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  const real Hn0_geom = Hn0_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
  return wvar(0, k + ks, j + js, i + is, n) * Hn0_geom;
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hn0(const real5d &uwvar, const real5d &wvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  real uw = compute_Hn0(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0, k + ks, j + js, i + is, n) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0, k + ks, j + js, i + is, n) += uw;
  }
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hn0(SArray<real, 1, 1> &uwvar, const real5d &wvar,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  real uw = compute_Hn0(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0) += uw;
  }
}

void YAKL_INLINE Hnm11bar_diagonal(SArray<real, 1, ndims> &Hnm11bardiag,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  for (int d = 0; d < ndims; d++) {
    Hnm11bardiag(d) = pgeom.get_area_10entity(d, k + ks, j + js, i + is, n) /
                      dgeom.get_area_nm11entity(d, k + ks, j + js, i + is, n);
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hnm11bar(SArray<real, 1, ndims> &u, const real5d &vvar,
                                  const Geometry<Straight> &pgeom,
                                  const Geometry<Twisted> &dgeom, int is,
                                  int js, int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;
  SArray<real, 1, ndims> Hnm11bargeom;

  Hnm11bar_diagonal(Hnm11bargeom, pgeom, dgeom, is, js, ks, i, j, k, n);

  for (int p = 0; p < ord - 1; p++) {
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        v(d, p) = vvar(d, k + ks, j + js, i + is + p - off, n);
      }
      if (d == 1) {
        v(d, p) = vvar(d, k + ks, j + js + p - off, i + is, n);
      }
    }
  }
  H1bar(u, v, Hnm11bargeom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hnm11bar(const real5d &uvar, const real5d &vvar,
                                  const Geometry<Straight> &pgeom,
                                  const Geometry<Twisted> &dgeom, int is,
                                  int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> u;
  compute_Hnm11bar<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k,
                                    n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      uvar(d, k + ks, j + js, i + is, n) = u(d);
    }
    if (addmode == ADD_MODE::ADD) {
      uvar(d, k + ks, j + js, i + is, n) += u(d);
    }
  }
}

real YAKL_INLINE Hn0bar_diagonal(const Geometry<Straight> &pgeom,
                                 const Geometry<Twisted> &dgeom, int is, int js,
                                 int ks, int i, int j, int k, int n) {
  return -pgeom.get_area_01entity(k + ks, j + js, i + is, n) /
         dgeom.get_area_n0entity(k + ks + 1, j + js, i + is, n);
}

real YAKL_INLINE compute_Hn0bar(const real5d &wvar,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  const real Hn0bargeom = Hn0bar_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
  return wvar(0, k + ks + 1, j + js, i + is, n) * Hn0bargeom;
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hn0bar(const real5d &uwvar, const real5d &wvar,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  real uw = compute_Hn0bar(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0, k + ks, j + js, i + is, n) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0, k + ks, j + js, i + is, n) += uw;
  }
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_Hn0bar(SArray<real, 1, 1> &uwvar, const real5d &wvar,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  real uw = compute_Hn0bar(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0) += uw;
  }
}

real YAKL_INLINE Hn1_diagonal(const Geometry<Straight> &pgeom,
                              const Geometry<Twisted> &dgeom, int is, int js,
                              int ks, int i, int j, int k, int n) {
  return dgeom.get_area_00entity(k + ks, j + js, i + is, n) /
         pgeom.get_area_n1entity(k + ks - 1, j + js, i + is, n);
}

// BROKEN FOR 2D+1D EXT
// JUST IN THE AREA CALCS...
template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_Hn1(SArray<real, 1, ndofs> &x0, const real5d &var,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real Hn1geom = Hn1_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int p = 0; p < hord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks - 1, j + js, i + is + p - hoff, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks - 1, j + js + p - hoff, i + is, n);
        }
      }
    }
  }
  H2<ndofs>(x0, x, Hn1geom);
  // EVENTUALLY BE MORE CLEVER IN THE VERTICAL HERE
  //  BUT THIS IS 2nd ORDER RIGHT NOW!
}

// Indexing is tricky, since we going from p11 to d00 and there is an "extended"
// boundary in d00 Since the i,j,k here are for d00, we must subtract 1 from k
// when indexing p11 quantities
//  ie the kth dual vertex corresponds with the k-1th primal cell
//  Also, should just compute over k=[1,..,ni-2] ie skip the first and last
//  vertices! These values should be set as boundary conditions... = 0 for no
//  flux I think!
template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_Hn1(const real5d &var0, const real5d &var,
                             const Geometry<Straight> &pgeom,
                             const Geometry<Twisted> &dgeom, int is, int js,
                             int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_Hn1<ndofs, hord, vord, hoff, voff>(x0, var, pgeom, dgeom, is, js, ks,
                                             i, j, k, n);

  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}
real YAKL_INLINE Hn1bar_diagonal(const Geometry<Straight> &pgeom,
                                 const Geometry<Twisted> &dgeom, int is, int js,
                                 int ks, int i, int j, int k, int n) {
  return pgeom.get_area_00entity(k + ks, j + js, i + is, n) /
         dgeom.get_area_n1entity(k + ks, j + js, i + is, n);
}

template <uint ndofs, uint vord, uint voff = vord / 2 - 1>
void YAKL_INLINE compute_Hn1bar(SArray<real, 1, ndofs> &x0, const real3d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int ks, int k,
                                int n) {
  const real Hn1bargeom = Hn1bar_diagonal(pgeom, dgeom, 0, 0, ks, 0, 0, k, n);
  for (int l = 0; l < ndofs; l++) {
    x0(l) = var(l, k + ks, n) * Hn1bargeom;
  }
}

template <uint ndofs, uint vord, uint voff = vord / 2 - 1, class F>
void YAKL_INLINE compute_Hn1bar(F f, SArray<real, 1, ndofs> &x0,
                                const real3d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int ks, int k,
                                int n) {
  const real Hn1bargeom = Hn1bar_diagonal(pgeom, dgeom, 0, 0, ks, 0, 0, k, n);
  for (int l = 0; l < ndofs; l++) {
    x0(l) = f(var, l, k + ks, n) * Hn1bargeom;
  }
}

// BROKEN FOR 2D+1D EXT
// JUST IN THE AREA FORM CALCS...
template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1, class F>
void YAKL_INLINE compute_Hn1bar(F f, SArray<real, 1, ndofs> &x0,
                                const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real Hn1bargeom = Hn1bar_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int p = 0; p < hord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        // not applying a transformation
        if constexpr (std::is_null_pointer_v<F>) {
          if (d == 0) {
            x(l, d, p) = var(l, k + ks, j + js, i + is + p - hoff, n);
          }
          if (d == 1) {
            x(l, d, p) = var(l, k + ks, j + js + p - hoff, i + is, n);
          }
        } else { // applying a transformation
          if (d == 0) {
            x(l, d, p) = f(var, l, k + ks, j + js, i + is + p - hoff, n);
          }
          if (d == 1) {
            x(l, d, p) = f(var, l, k + ks, j + js + p - hoff, i + is, n);
          }
        }
      }
    }
  }
  H2bar<ndofs>(x0, x, Hn1bargeom);
  // EVENTUALLY BE MORE CLEVER IN THE VERTICAL HERE
  //  BUT THIS IS 2nd ORDER RIGHT NOW!
}
template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_Hn1bar(SArray<real, 1, ndofs> &x0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  compute_Hn1bar<ndofs, hord, vord>(nullptr, x0, var, pgeom, dgeom, is, js, ks,
                                    i, j, k, n);
}

// Indexing here is fine, since we going from d11 to p00 and there is no
// "extended" boundary in p00
template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1, class F>
void YAKL_INLINE compute_Hn1bar(F f, const real5d &var0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_Hn1bar<ndofs, hord, vord, hoff, voff>(f, x0, var, pgeom, dgeom, is,
                                                js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}
template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_Hn1bar(const real5d &var0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  compute_Hn1bar<ndofs, hord, vord, addmode, hoff, voff>(
      nullptr, var0, var, pgeom, dgeom, is, js, ks, i, j, k, n);
}

template <uint ord, int off = ord / 2 - 1>
real YAKL_INLINE fourier_Hn1bar(const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n, int nx,
                                int ny, int nz) {
  SArray<real, 2, ndims, off + GPU_PAD> shift;

  // assuming these are constant
  const real Hn1bargeom = Hn1bar_diagonal(pgeom, dgeom, is, js, ks, i, j, k, n);

  for (int p = 0; p < off; p++) {
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        shift(d, p) = (2 * pi * i * (p - off)) / nx;
      }
      if (d == 1) {
        shift(d, p) = (2 * pi * j * (p - off)) / ny;
      }
    }
  }
  return H2barhat(Hn1bargeom, shift);
}
} // namespace pamc
