
#pragma once

#include "common.h"
#include "geometry.h"

void YAKL_INLINE H1(SArray<real, 1, ndims> &var,
                    SArray<real, 2, ndims, 1> const &velocity,
                    SArray<real, 1, ndims> const &H1geom) {

  for (int d = 0; d < ndims; d++) {
    var(d) = velocity(d, 0);
    var(d) *= H1geom(d);
  }
}

void YAKL_INLINE H1(SArray<real, 1, ndims> &var,
                    SArray<real, 2, ndims, 3> const &velocity,
                    SArray<real, 1, ndims> const &H1geom) {

  for (int d = 0; d < ndims; d++) {
    var(d) = -1.0_fp / 24.0_fp * velocity(d, 0) +
             26.0_fp / 24.0_fp * velocity(d, 1) -
             1.0_fp / 24.0_fp * velocity(d, 2);
    var(d) *= H1geom(d);
  }
}

void YAKL_INLINE H1(SArray<real, 1, ndims> &var,
                    SArray<real, 2, ndims, 5> const &velocity,
                    SArray<real, 1, ndims> const &H1geom) {

  for (int d = 0; d < ndims; d++) {
    var(d) = 9.0_fp / 1920.0_fp * velocity(d, 0) -
             116.0_fp / 1920.0_fp * velocity(d, 1) +
             2134.0_fp / 1920.0_fp * velocity(d, 2) -
             116.0_fp / 1920.0_fp * velocity(d, 3) +
             9.0_fp / 1920.0_fp * velocity(d, 4);

    var(d) *= H1geom(d);
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1(SArray<real, 1, ndims> &u, const real5d &vvar,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;
  SArray<real, 1, ndims> H1geom;
  for (int d = 0; d < ndims; d++) {
    H1geom(d) = dgeom.get_area_lform(ndims - 1, d, k + ks, j + js, i + is, n) /
                pgeom.get_area_lform(1, d, k + ks, j + js, i + is, n);
  }

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
  H1(u, v, H1geom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1(const real5d &uvar, const real5d &vvar,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> u;
  compute_H1<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    for (int d = 0; d < ndims; d++) {
      uvar(d, k + ks, j + js, i + is, n) = u(d);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int d = 0; d < ndims; d++) {
      uvar(d, k + ks, j + js, i + is, n) += u(d);
    }
  }
}

void YAKL_INLINE H1hat(SArray<real, 1, ndims> &u,
                       SArray<real, 1, ndims> const &H1geom,
                       SArray<real, 2, ndims, 1> const &shift) {
  for (int d = 0; d < ndims; d++) {
    u(d) = H1geom(d);
  }
}

void YAKL_INLINE H1hat(SArray<real, 1, ndims> &u,
                       SArray<real, 1, ndims> const &H1geom,
                       SArray<real, 2, ndims, 2> const &shift) {
  for (int d = 0; d < ndims; d++) {
    u(d) = -2.0_fp / 24.0_fp * cos(shift(d, 0)) + 26.0_fp / 24.0_fp;
    u(d) *= H1geom(d);
  }
}

void YAKL_INLINE H1hat(SArray<real, 1, ndims> &u,
                       SArray<real, 1, ndims> const &H1geom,
                       SArray<real, 2, ndims, 3> const &shift) {
  for (int d = 0; d < ndims; d++) {
    u(d) = 18.0_fp / 1920.0_fp * cos(shift(d, 0)) -
           232.0_fp / 1920.0_fp * cos(shift(d, 1)) + 2134.0_fp / 1920.0_fp;
    u(d) *= H1geom(d);
  }
}

template <uint ord, int off = ord / 2 - 1>
void YAKL_INLINE fourier_H1(SArray<real, 1, ndims> &u,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n, int nx, int ny,
                            int nz) {
  // Adding 1 to off because GPUs don't like zero size structs
  SArray<real, 2, ndims, off + 1> shift;
  SArray<real, 1, ndims> H1geom;
  for (int d = 0; d < ndims; d++) {
    H1geom(d) = dgeom.get_area_lform(ndims - 1, d, k + ks, j + js, i + is, n) /
                pgeom.get_area_lform(1, d, k + ks, j + js, i + is, n);
  }

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
  H1hat(u, H1geom, shift);
}

// Note the indexing here, this is key
real YAKL_INLINE compute_H1_vert(const real5d &wvar,
                                 const Geometry<Straight> &pgeom,
                                 const Geometry<Twisted> &dgeom, int is, int js,
                                 int ks, int i, int j, int k, int n) {
  // THIS IS 2ND ORDER AT BEST...
  return wvar(0, k + ks - 1, j + js, i + is, n) *
         dgeom.get_area_10entity(k + ks, j + js, i + is, n) /
         pgeom.get_area_01entity(k + ks - 1, j + js, i + is, n);
}

real YAKL_INLINE H1_vert_coeff(const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k) {
  // THIS IS 2ND ORDER AT BEST...
  return dgeom.get_area_10entity(k + ks, j + js, i + is, n) /
         pgeom.get_area_01entity(k + ks - 1, j + js, i + is, n);
}

// Indexing issues since we go from p01 to d10, and d10 has an "extended
// boundary" ie boundary vert edges Since we index over d10, need to subtract 1
// from k when indexing p01 ie the kth edge flux corresponds with the k-1th edge
// velocity Also should be called with k=[1,...,ni-2] ie skip the first and last
// fluxes, which are set diagnostically (=0 for no-flux bcs)
template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1_vert(const real5d &uwvar, const real5d &wvar,
                                 const Geometry<Straight> &pgeom,
                                 const Geometry<Twisted> &dgeom, int is, int js,
                                 int ks, int i, int j, int k, int n) {
  real uw = compute_H1_vert(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0, k + ks, j + js, i + is, n) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0, k + ks, j + js, i + is, n) += uw;
  }
}
template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1_vert(SArray<real, 1, 1> &uwvar, const real5d &wvar,
                                 const Geometry<Straight> &pgeom,
                                 const Geometry<Twisted> &dgeom, int is, int js,
                                 int ks, int i, int j, int k, int n) {
  real uw = compute_H1_vert(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0) += uw;
  }
}

// BROKEN FOR 2D+1D EXT
// MAINLY IN THE AREA CALCS...
template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1_ext(SArray<real, 1, ndims> &u, const real5d &vvar,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;
  SArray<real, 1, ndims> H1geom;

  for (int d = 0; d < ndims; d++) {
    H1geom(d) = dgeom.get_area_01entity(k + ks, j + js, i + is, n) /
                pgeom.get_area_10entity(k + ks, j + js, i + is, n);
  }

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
  H1(u, v, H1geom);
}

// No indexing issues since we go from p10 to d01, and d01 has no "extended
// boundary"
template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1_ext(const real5d &uvar, const real5d &vvar,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> u;
  compute_H1_ext<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k,
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

template <uint ord, int off = ord / 2 - 1>
void YAKL_INLINE fourier_H1_ext(SArray<real, 1, ndims> &u,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n, int nx,
                                int ny, int nz) {

  // Adding 1 to off because GPUs don't like zero size structs
  SArray<real, 2, ndims, off + 1> shift;
  SArray<real, 1, ndims> H1geom;
  for (int d = 0; d < ndims; d++) {
    H1geom(d) = dgeom.get_area_01entity(k + ks, j + js, i + is, n) /
                pgeom.get_area_10entity(k + ks, j + js, i + is, n);
  }

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
  H1hat(u, H1geom, shift);
}

template <uint ndofs>
void YAKL_INLINE H2bar(SArray<real, 1, ndofs> &var,
                       SArray<real, 3, ndofs, ndims, 1> const &dens,
                       real H2bargeom) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 0);
    var(l) *= H2bargeom;
  }
}

template <uint ndofs>
void YAKL_INLINE H2bar(SArray<real, 1, ndofs> &var,
                       SArray<real, 3, ndofs, ndims, 3> const &dens,
                       real H2bargeom) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 1);
    for (int d = 0; d < ndims; d++) {
      var(l) += -1.0_fp / 24.0_fp * dens(l, d, 0) +
                2.0_fp / 24.0_fp * dens(l, d, 1) -
                1.0_fp / 24.0_fp * dens(l, d, 2);
    }
    var(l) *= H2bargeom;
  }
}

template <uint ndofs>
void YAKL_INLINE H2bar(SArray<real, 1, ndofs> &var,
                       SArray<real, 3, ndofs, ndims, 5> const &dens,
                       real H2bargeom) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 2);
    for (int d = 0; d < ndims; d++) {
      var(l) += 9.0_fp / 1920.0_fp * dens(l, d, 0) -
                116.0_fp / 1920.0_fp * dens(l, d, 1) +
                214.0_fp / 1920.0_fp * dens(l, d, 2) -
                116.0_fp / 1920.0_fp * dens(l, d, 3) +
                9.0_fp / 1920.0_fp * dens(l, d, 4);
    }
    var(l) *= H2bargeom;
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H2bar(SArray<real, 1, ndofs> &x0, const real5d &var,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H2bargeom = pgeom.get_area_lform(0, 0, k + ks, j + js, i + is, n) /
                         dgeom.get_area_lform(ndims, 0, k + ks, j + js, i + is, n);

  for (int p = 0; p < ord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks, j + js, i + is + p - off, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks, j + js + p - off, i + is, n);
        }
      }
    }
  }
  H2bar<ndofs>(x0, x, H2bargeom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H2bar(const real5d &var0, const real5d &var,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H2bar<ndofs, ord, off>(x0, var, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}

real YAKL_INLINE H2barhat(real H2bargeom,
                          const SArray<real, 2, ndims, 1> &shift) {
  return H2bargeom;
}
real YAKL_INLINE H2barhat(real &H2bargeom,
                          const SArray<real, 2, ndims, 2> &shift) {
  real res = 1.0_fp;
  for (int d = 0; d < ndims; d++) {
    res += -2.0_fp / 24.0_fp * cos(shift(d, 0)) + 2.0_fp / 24.0_fp;
  }
  return H2bargeom * res;
}
real YAKL_INLINE H2barhat(real H2bargeom,
                          const SArray<real, 2, ndims, 3> &shift) {
  real res = 1.0_fp;
  for (int d = 0; d < ndims; d++) {
    res += 18.0_fp / 1920.0_fp * cos(shift(d, 0)) -
           232.0_fp / 1920.0_fp * cos(shift(d, 1)) + 214.0_fp / 1920.0_fp;
  }
  return H2bargeom * res;
}

template <uint ord, int off = ord / 2 - 1>
real YAKL_INLINE fourier_H2bar(const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n, int nx,
                               int ny, int nz) {
  // Adding 1 to off because GPUs don't like zero size structs
  SArray<real, 2, ndims, off + 1> shift;

  // assuming these are constant
  const real H2bargeom = pgeom.get_area_lform(0, 0, k + ks, j + js, i + is, n) /
                         dgeom.get_area_lform(ndims, 0, k + ks, j + js, i + is, n);

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
  return H2barhat(H2bargeom, shift);
}

template <uint ndofs, uint vord, uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H2bar_ext(SArray<real, 1, ndofs> &x0,
                                   const real3d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int ks,
                                   int k, int n) {
  real H2bargeom = pgeom.get_area_00entity(k + ks, 0, 0, n) /
                   dgeom.get_area_11entity(k + ks, 0, 0, n);
  for (int l = 0; l < ndofs; l++) {
    x0(l) = var(l, k, n) * H2bargeom;
  }
}

// BROKEN FOR 2D+1D EXT
// JUST IN THE AREA FORM CALCS...
template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1, class F>
void YAKL_INLINE compute_H2bar_ext(F f, SArray<real, 1, ndofs> &x0,
                                   const real5d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real H2bargeom = pgeom.get_area_00entity(k + ks, j + js, i + is, n) /
                         dgeom.get_area_11entity(k + ks, j + js, i + is, n);
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
  H2bar<ndofs>(x0, x, H2bargeom);
  // EVENTUALLY BE MORE CLEVER IN THE VERTICAL HERE
  //  BUT THIS IS 2nd ORDER RIGHT NOW!
}
template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H2bar_ext(SArray<real, 1, ndofs> &x0,
                                   const real5d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  compute_H2bar_ext<ndofs, hord, vord>(nullptr, x0, var, pgeom, dgeom, is, js,
                                       ks, i, j, k, n);
}

// Indexing here is fine, since we going from d11 to p00 and there is no
// "extended" boundary in p00
template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1, class F>
void YAKL_INLINE compute_H2bar_ext(F f, const real5d &var0, const real5d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H2bar_ext<ndofs, hord, vord, hoff, voff>(f, x0, var, pgeom, dgeom, is,
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
void YAKL_INLINE compute_H2bar_ext(const real5d &var0, const real5d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  compute_H2bar_ext<ndofs, hord, vord, addmode, hoff, voff>(
      nullptr, var0, var, pgeom, dgeom, is, js, ks, i, j, k, n);
}

template <uint ord, int off = ord / 2 - 1>
real YAKL_INLINE fourier_H2bar_ext(const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n,
                                   int nx, int ny, int nz) {
  // Adding 1 to off because GPUs don't like zero size structs
  SArray<real, 2, ndims, off + 1> shift;

  // assuming these are constant
  const real H2bargeom =
      pgeom.get_area_00entity(k + ks, j + js, i + is + 0 - off, n) /
      dgeom.get_area_11entity(k + ks, j + js, i + is + 0 - off, n);

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
  return H2barhat(H2bargeom, shift);
}

template <uint ndofs>
void YAKL_INLINE H2(SArray<real, 1, ndofs> &var,
                    SArray<real, 3, ndofs, ndims, 1> const &dens, real H2geom) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 0);
    var(l) *= H2geom;
  }
}

template <uint ndofs>
void YAKL_INLINE H2(SArray<real, 1, ndofs> &var,
                    SArray<real, 3, ndofs, ndims, 3> const &dens, real H2geom) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 1);
    for (int d = 0; d < ndims; d++) {
      var(l) += -1.0_fp / 24.0_fp * dens(l, d, 0) +
                2.0_fp / 24.0_fp * dens(l, d, 1) -
                1.0_fp / 24.0_fp * dens(l, d, 2);
    }
    var(l) *= H2geom;
  }
}

template <uint ndofs>
void YAKL_INLINE H2(SArray<real, 1, ndofs> &var,
                    SArray<real, 3, ndofs, ndims, 5> const &dens, real H2geom) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 2);
    for (int d = 0; d < ndims; d++) {
      var(l) += 9.0_fp / 1920.0_fp * dens(l, d, 0) -
                116.0_fp / 1920.0_fp * dens(l, d, 1) +
                214.0_fp / 1920.0_fp * dens(l, d, 2) -
                116.0_fp / 1920.0_fp * dens(l, d, 3) +
                9.0_fp / 1920.0_fp * dens(l, d, 4);
    }
    var(l) *= H2geom;
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H2(SArray<real, 1, ndofs> &x0, const real5d &var,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H2geom = dgeom.get_area_lform(0, 0, k + ks, j + js, i + is, n) /
                      pgeom.get_area_lform(ndims, 0, k + ks, j + js, i + is, n);
  for (int p = 0; p < ord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks, j + js, i + is + p - off, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks, j + js + p - off, i + is, n);
        }
      }
    }
  }
  H2<ndofs>(x0, x, H2geom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H2(const real5d &var0, const real5d &var,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H2<ndofs, ord, off>(x0, var, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}

// BROKEN FOR 2D+1D EXT
// JUST IN THE AREA CALCS...
template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H2_ext(SArray<real, 1, ndofs> &x0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real H2geom = dgeom.get_area_00entity(k + ks, j + js, i + is, n) /
                      pgeom.get_area_11entity(k + ks - 1, j + js, i + is, n);
  for (int p = 0; p < hord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks - 1, j + js, i + is + p - hoff, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks + 1, j + js + p - hoff, i + is, n);
        }
      }
    }
  }
  H2<ndofs>(x0, x, H2geom);
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
void YAKL_INLINE compute_H2_ext(const real5d &var0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H2_ext<ndofs, hord, vord, hoff, voff>(x0, var, pgeom, dgeom, is, js,
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

template <uint ndofs>
void YAKL_INLINE H0(SArray<real, 1, ndofs> &var,
                    SArray<real, 3, ndofs, ndims, 1> const &dens, real H0geom) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 0);
    var(l) *= H0geom;
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H0(SArray<real, 1, ndofs> &x0, const real5d &var,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H0geom = dgeom.get_area_lform(2, 0, k + ks, j + js, i + is, n) /
                      pgeom.get_area_lform(0, 0, k + ks, j + js, i + is, n);
  for (int p = 0; p < ord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks, j + js, i + is + p - off, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks, j + js + p - off, i + is, n);
        }
      }
    }
  }
  H0<ndofs>(x0, x, H0geom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H0(const real5d &var0, const real5d &var,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H0<ndofs, ord, off>(x0, var, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}

template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H0_ext(SArray<real, 1, ndofs> &x0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real H0geom = dgeom.get_area_11entity(k + ks, j + js, i + is, n) /
                      pgeom.get_area_00entity(k + ks, j + js, i + is, n);
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
  H0<ndofs>(x0, x, H0geom);
}

template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H0_ext(const real5d &var0, const real5d &var,
                                const Geometry<Straight> &pgeom,
                                const Geometry<Twisted> &dgeom, int is, int js,
                                int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H0_ext<ndofs, hord, vord, hoff, voff>(x0, var, pgeom, dgeom, is, js,
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

template <uint ndofs>
void YAKL_INLINE H0bar(SArray<real, 1, ndofs> &var,
                       SArray<real, 3, ndofs, ndims, 1> const &dens,
                       real H0bargeom) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 0);
    var(l) *= H0bargeom;
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H0bar(SArray<real, 1, ndofs> &x0, const real5d &var,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H0bargeom = dgeom.get_area_lform(2, 0, k + ks, j + js, i + is, n) /
                         pgeom.get_area_lform(0, 0, k + ks, j + js, i + is, n);
  for (int p = 0; p < ord - 1; p++) {
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          x(l, d, p) = var(l, k + ks, j + js, i + is + p - off, n);
        }
        if (d == 1) {
          x(l, d, p) = var(l, k + ks, j + js + p - off, i + is, n);
        }
      }
    }
  }
  H0bar<ndofs>(x0, x, H0bargeom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H0bar(const real5d &var0, const real5d &var,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H0bar<ndofs, ord, off>(x0, var, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      var0(l, k + ks, j + js, i + is, n) = x0(l);
    }
    if (addmode == ADD_MODE::ADD) {
      var0(l, k + ks, j + js, i + is, n) += x0(l);
    }
  }
}

template <uint ndofs, uint hord, uint vord, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H0bar_ext(SArray<real, 1, ndofs> &x0,
                                   const real5d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, hord - 1> x;
  const real H0bargeom = pgeom.get_area_11entity(k + ks, j + js, i + is, n) /
                         dgeom.get_area_00entity(k + ks, j + js, i + is, n);
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
  H0bar<ndofs>(x0, x, H0bargeom);
}

template <uint ndofs, uint hord, uint vord,
          ADD_MODE addmode = ADD_MODE::REPLACE, uint hoff = hord / 2 - 1,
          uint voff = vord / 2 - 1>
void YAKL_INLINE compute_H0bar_ext(const real5d &var0, const real5d &var,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> x0;
  compute_H0bar_ext<ndofs, hord, vord, hoff, voff>(x0, var, pgeom, dgeom, is,
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

void YAKL_INLINE H1bar(SArray<real, 1, ndims> &var,
                       SArray<real, 2, ndims, 1> const &velocity,
                       SArray<real, 1, ndims> const &H1bargeom) {

  for (int d = 0; d < ndims; d++) {
    var(d) = -velocity(d, 0);
    var(d) *= H1bargeom(d);
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar(SArray<real, 1, ndims> &u, const real5d &vvar,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;
  SArray<real, 1, ndims> H1bargeom;
  for (int d = 0; d < ndims; d++) {
    H1bargeom(d) = pgeom.get_area_lform(1, d, k + ks, j + js, i + is, n) /
                   dgeom.get_area_lform(ndims - 1, d, k + ks, j + js, i + is, n);
  }

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
  H1bar(u, v, H1bargeom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar(const real5d &uvar, const real5d &vvar,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> u;
  compute_H1bar<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    for (int d = 0; d < ndims; d++) {
      uvar(d, k + ks, j + js, i + is, n) = u(d);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int d = 0; d < ndims; d++) {
      uvar(d, k + ks, j + js, i + is, n) += u(d);
    }
  }
}

real YAKL_INLINE compute_H1bar_vert(const real5d &wvar,
                                    const Geometry<Straight> &pgeom,
                                    const Geometry<Twisted> &dgeom, int is,
                                    int js, int ks, int i, int j, int k,
                                    int n) {
  return -wvar(0, k + ks + 1, j + js, i + is, n) *
         pgeom.get_area_01entity(k + ks, j + js, i + is, n) /
         dgeom.get_area_10entity(k + ks + 1, j + js, i + is, n);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar_vert(const real5d &uwvar, const real5d &wvar,
                                    const Geometry<Straight> &pgeom,
                                    const Geometry<Twisted> &dgeom, int is,
                                    int js, int ks, int i, int j, int k,
                                    int n) {
  real uw = compute_H1bar_vert(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0, k + ks, j + js, i + is, n) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0, k + ks, j + js, i + is, n) += uw;
  }
}
template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar_vert(SArray<real, 1, 1> &uwvar,
                                    const real5d &wvar,
                                    const Geometry<Straight> &pgeom,
                                    const Geometry<Twisted> &dgeom, int is,
                                    int js, int ks, int i, int j, int k,
                                    int n) {
  real uw = compute_H1bar_vert(wvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    uwvar(0) = uw;
  }
  if (addmode == ADD_MODE::ADD) {
    uwvar(0) += uw;
  }
}

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar_ext(SArray<real, 1, ndims> &u,
                                   const real5d &vvar,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;
  SArray<real, 1, ndims> H1bargeom;

  for (int d = 0; d < ndims; d++) {
    H1bargeom(d) = pgeom.get_area_10entity(k + ks, j + js, i + is, n) /
                   dgeom.get_area_01entity(k + ks, j + js, i + is, n);
  }

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
  H1bar(u, v, H1bargeom);
}

template <uint ndofs, uint ord, ADD_MODE addmode = ADD_MODE::REPLACE,
          uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar_ext(const real5d &uvar, const real5d &vvar,
                                   const Geometry<Straight> &pgeom,
                                   const Geometry<Twisted> &dgeom, int is,
                                   int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> u;
  compute_H1bar_ext<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k,
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
