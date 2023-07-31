#pragma once

#include "common.h"
#include "geometry.h"

namespace pamc {

template <uint ndofs>
void YAKL_INLINE H0(SArray<real, 1, ndofs> &var,
                    SArray<real, 3, ndofs, ndims, 1> const &dens, real H0geom) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = dens(l, 0, 0);
    var(l) *= H0geom;
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

void YAKL_INLINE H1bar(SArray<real, 1, ndims> &var,
                       SArray<real, 2, ndims, 1> const &velocity,
                       SArray<real, 1, ndims> const &H1bargeom) {

  for (int d = 0; d < ndims; d++) {
    var(d) = -velocity(d, 0);
    var(d) *= H1bargeom(d);
  }
}

void YAKL_INLINE H1hat(SArray<real, 1, ndims> &u,
                       SArray<real, 1, ndims> const &H1geom,
                       SArray<real, 2, ndims, GPU_PAD> const &shift) {
  for (int d = 0; d < ndims; d++) {
    u(d) = H1geom(d);
  }
}

void YAKL_INLINE H1hat(SArray<real, 1, ndims> &u,
                       SArray<real, 1, ndims> const &H1geom,
                       SArray<real, 2, ndims, 1 + GPU_PAD> const &shift) {
  for (int d = 0; d < ndims; d++) {
    u(d) = -2.0_fp / 24.0_fp * cos(shift(d, 0)) + 26.0_fp / 24.0_fp;
    u(d) *= H1geom(d);
  }
}

void YAKL_INLINE H1hat(SArray<real, 1, ndims> &u,
                       SArray<real, 1, ndims> const &H1geom,
                       SArray<real, 2, ndims, 2 + GPU_PAD> const &shift) {
  for (int d = 0; d < ndims; d++) {
    u(d) = 18.0_fp / 1920.0_fp * cos(shift(d, 0)) -
           232.0_fp / 1920.0_fp * cos(shift(d, 1)) + 2134.0_fp / 1920.0_fp;
    u(d) *= H1geom(d);
  }
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

real YAKL_INLINE H2barhat(real H2bargeom,
                          const SArray<real, 2, ndims, GPU_PAD> &shift) {
  return H2bargeom;
}

real YAKL_INLINE H2barhat(real &H2bargeom,
                          const SArray<real, 2, ndims, 1 + GPU_PAD> &shift) {
  real res = 1.0_fp;
  for (int d = 0; d < ndims; d++) {
    res += -2.0_fp / 24.0_fp * cos(shift(d, 0)) + 2.0_fp / 24.0_fp;
  }
  return H2bargeom * res;
}

real YAKL_INLINE H2barhat(real H2bargeom,
                          const SArray<real, 2, ndims, 2 + GPU_PAD> &shift) {
  real res = 1.0_fp;
  for (int d = 0; d < ndims; d++) {
    res += 18.0_fp / 1920.0_fp * cos(shift(d, 0)) -
           232.0_fp / 1920.0_fp * cos(shift(d, 1)) + 214.0_fp / 1920.0_fp;
  }
  return H2bargeom * res;
}
} // namespace pamc

#include "hodge_star_extruded.h"
#include "hodge_star_layer.h"
