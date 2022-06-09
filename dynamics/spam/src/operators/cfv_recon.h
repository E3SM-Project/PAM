#pragma once

#include "common.h"

template <uint ndofs, uint nd>
void YAKL_INLINE cfv(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                     SArray<real, 3, ndofs, nd, 1> const &dens) {
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      edgerecon(l, d, 0) = dens(l, d, 0);
      edgerecon(l, d, 1) = dens(l, d, 0);
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE cfv(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                     SArray<real, 3, ndofs, nd, 3> const &dens) {
  real er;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      er = (8.0_fp / 6.0_fp) * dens(l, d, 1) -
           (1.0_fp / 6.0_fp) * (dens(l, d, 0) + dens(l, d, 2));
      edgerecon(l, d, 0) = er;
      edgerecon(l, d, 1) = er;
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE cfv(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                     SArray<real, 3, ndofs, nd, 5> const &dens) {
  real er;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      er = (46.0_fp / 30.0_fp) * dens(l, d, 2) -
           (9.0_fp / 30.0_fp) * (dens(l, d, 1) + dens(l, d, 3)) +
           (1.0_fp / 30.0_fp) * (dens(l, d, 0) + dens(l, d, 4));
      edgerecon(l, d, 0) = er;
      edgerecon(l, d, 1) = er;
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE cfv(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                     SArray<real, 3, ndofs, nd, 7> const &dens) {
  real er;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      er = (704.0_fp / 420.0_fp) * dens(l, d, 3) -
           (171.0_fp / 420.0_fp) * (dens(l, d, 2) + dens(l, d, 4)) +
           (32.0_fp / 420.0_fp) * (dens(l, d, 1) + dens(l, d, 5)) -
           (3.0_fp / 420.0_fp) * (dens(l, d, 0) + dens(l, d, 6));
      edgerecon(l, d, 0) = er;
      edgerecon(l, d, 1) = er;
    }
  }
}

template <uint ndofs, uint nd>
void YAKL_INLINE cfv(SArray<real, 3, ndofs, nd, 2> &edgerecon,
                     SArray<real, 3, ndofs, nd, 9> const &dens) {
  real er;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < nd; d++) {
      er = (2252.0_fp / 1260.0_fp) * dens(l, d, 4) -
           (625.0_fp / 1260.0_fp) * (dens(l, d, 3) + dens(l, d, 5)) +
           (152.0_fp / 1260.0_fp) * (dens(l, d, 2) + dens(l, d, 6)) -
           (25.0_fp / 1260.0_fp) * (dens(l, d, 1) + dens(l, d, 7)) +
           (2.0_fp / 1260.0_fp) * (dens(l, d, 0) + dens(l, d, 8));
      edgerecon(l, d, 0) = er;
      edgerecon(l, d, 1) = er;
    }
  }
}