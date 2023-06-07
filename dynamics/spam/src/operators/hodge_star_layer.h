#pragma once

#include "common.h"
#include "geometry.h"
#include "hodge_star.h"

namespace pamc {

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H0(SArray<real, 1, ndofs> &x0, const real5d &var,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H0geom =
      dgeom.get_area_lform<ndims>(0, k + ks, j + js, i + is, n) /
      pgeom.get_area_lform<0>(0, k + ks, j + js, i + is, n);
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

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H0bar(SArray<real, 1, ndofs> &x0, const real5d &var,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H0bargeom =
      dgeom.get_area_lform<ndims>(0, k + ks, j + js, i + is, n) /
      pgeom.get_area_lform<0>(0, k + ks, j + js, i + is, n);
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

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1(SArray<real, 2, ndofs, ndims> &u,
                            const real5d &vvar, const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> v;
  SArray<real, 1, ndims> H1geom;
  for (int d = 0; d < ndims; d++) {
    H1geom(d) = dgeom.get_area_lform<ndims - 1>(d, k + ks, j + js, i + is, n) /
                pgeom.get_area_lform<1>(d, k + ks, j + js, i + is, n);
  }

  for (int l = 0; l < ndofs; l++) {
    for (int p = 0; p < ord - 1; p++) {
      for (int d = 0; d < ndims; d++) {
        if (d == 0) {
          v(l, d, p) = vvar(d, k + ks, j + js, i + is + p - off, n);
        }
        if (d == 1) {
          v(l, d, p) = vvar(d, k + ks, j + js + p - off, i + is, n);
        }
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
  SArray<real, 2, ndofs, ndims> u;
  compute_H1<ndofs, ord, off>(u, vvar, pgeom, dgeom, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      for (int d = 0; d < ndims; d++) {
        uvar(l * ndims + d, k + ks, j + js, i + is, n) = u(l, d);
      }
    }
    if (addmode == ADD_MODE::ADD) {
      for (int d = 0; d < ndims; d++) {
        uvar(l * ndims + d, k + ks, j + js, i + is, n) += u(l, d);
      }
    }
  }
}

template <uint ord, int off = ord / 2 - 1>
void YAKL_INLINE fourier_H1(SArray<real, 1, ndims> &u,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n, int nx, int ny,
                            int nz) {
  SArray<real, 2, ndims, off + GPU_PAD> shift;
  SArray<real, 1, ndims> H1geom;
  for (int d = 0; d < ndims; d++) {
    H1geom(d) = dgeom.get_area_lform<ndims - 1>(d, k + ks, j + js, i + is, n) /
                pgeom.get_area_lform<1>(d, k + ks, j + js, i + is, n);
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

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H1bar(SArray<real, 1, ndims> &u, const real5d &vvar,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, ord - 1> v;
  SArray<real, 1, ndims> H1bargeom;
  for (int d = 0; d < ndims; d++) {
    H1bargeom(d) =
        pgeom.get_area_lform<1>(d, k + ks, j + js, i + is, n) /
        dgeom.get_area_lform<ndims - 1>(d, k + ks, j + js, i + is, n);
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

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H2(SArray<real, 1, ndofs> &x0, const real5d &var,
                            const Geometry<Straight> &pgeom,
                            const Geometry<Twisted> &dgeom, int is, int js,
                            int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H2geom = dgeom.get_area_lform<0>(0, k + ks, j + js, i + is, n) /
                      pgeom.get_area_lform<ndims>(0, k + ks, j + js, i + is, n);
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

template <uint ndofs, uint ord, uint off = ord / 2 - 1>
void YAKL_INLINE compute_H2bar(SArray<real, 1, ndofs> &x0, const real5d &var,
                               const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, ord - 1> x;
  const real H2bargeom =
      pgeom.get_area_lform<0>(0, k + ks, j + js, i + is, n) /
      dgeom.get_area_lform<ndims>(0, k + ks, j + js, i + is, n);

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

template <uint ord, int off = ord / 2 - 1>
real YAKL_INLINE fourier_H2bar(const Geometry<Straight> &pgeom,
                               const Geometry<Twisted> &dgeom, int is, int js,
                               int ks, int i, int j, int k, int n, int nx,
                               int ny, int nz) {
  SArray<real, 2, ndims, off + GPU_PAD> shift;

  // assuming these are constant
  const real H2bargeom =
      pgeom.get_area_lform<0>(0, k + ks, j + js, i + is, n) /
      dgeom.get_area_lform<ndims>(0, k + ks, j + js, i + is, n);

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
} // namespace pamc
