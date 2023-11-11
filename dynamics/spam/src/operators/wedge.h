#pragma once

#include "common.h"

namespace pamc {

// Q
template <index_t ndofs>
void YAKL_INLINE Q2D(SArray<real, 2, ndofs, 2> &vel,
                     SArray<real, 3, ndofs, 2, 5> const &recon,
                     SArray<real, 2, 2, 4> const &flux) {

  for (int l = 0; l < ndofs; l++) {
    // x-dir
    vel(l, 0) =
        -0.125_fp * (flux(0, 0) * recon(l, 0, 0) + flux(0, 1) * recon(l, 0, 1) +
                     flux(0, 2) * recon(l, 0, 2) + flux(0, 3) * recon(l, 0, 3) +
                     flux(0, 0) * recon(l, 0, 4) + flux(0, 1) * recon(l, 0, 4) +
                     flux(0, 2) * recon(l, 0, 4) + flux(0, 3) * recon(l, 0, 4));

    // y-dir
    vel(l, 1) =
        0.125_fp * (flux(1, 0) * recon(l, 1, 0) + flux(1, 1) * recon(l, 1, 1) +
                    flux(1, 2) * recon(l, 1, 2) + flux(1, 3) * recon(l, 1, 3) +
                    flux(1, 0) * recon(l, 1, 4) + flux(1, 1) * recon(l, 1, 4) +
                    flux(1, 2) * recon(l, 1, 4) + flux(1, 3) * recon(l, 1, 4));
  }
}

template <index_t ndofs>
void YAKL_INLINE Q2D_nonEC(SArray<real, 2, ndofs, 2> &vel,
                           SArray<real, 2, ndofs, 2> const &recon,
                           SArray<real, 2, 2, 4> const &flux) {
  for (int l = 0; l < ndofs; l++) {
    // x-dir
    vel(l, 0) = -0.25_fp * recon(l, 0) *
                (flux(0, 0) + flux(0, 1) + flux(0, 2) + flux(0, 3));
    // y-dir
    vel(l, 1) = 0.25_fp * recon(l, 1) *
                (flux(1, 0) + flux(1, 1) + flux(1, 2) + flux(1, 3));
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Q_EC(const real5d &qflux, const real5d &reconvar,
                              const real5d &Uvar, int is, int js, int ks, int i,
                              int j, int k, int n) {

  SArray<real, 2, ndofs, 2> vel;
  SArray<real, 2, 2, 4> flux;
  flux(0, 0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(0, 1) = Uvar(1, k + ks, j + js, i + is - 1, n);
  flux(0, 2) = Uvar(1, k + ks, j + js + 1, i + is, n);
  flux(0, 3) = Uvar(1, k + ks, j + js + 1, i + is - 1, n);

  flux(1, 0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1, 1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  flux(1, 2) = Uvar(0, k + ks, j + js - 1, i + is, n);
  flux(1, 3) = Uvar(0, k + ks, j + js - 1, i + is + 1, n);

  SArray<real, 3, ndofs, 2, 5> recon;
  for (int l = 0; l < ndofs; l++) {
    recon(l, 0, 0) = reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    recon(l, 0, 1) = reconvar(l + 1 * ndofs, k + ks, j + js, i + is - 1, n);
    recon(l, 0, 2) = reconvar(l + 1 * ndofs, k + ks, j + js + 1, i + is, n);
    recon(l, 0, 3) = reconvar(l + 1 * ndofs, k + ks, j + js + 1, i + is - 1, n);
    recon(l, 0, 4) = reconvar(l + 0 * ndofs, k + ks, j + js, i + is, n);

    recon(l, 1, 0) = reconvar(l + 0 * ndofs, k + ks, j + js, i + is, n);
    recon(l, 1, 1) = reconvar(l + 0 * ndofs, k + ks, j + js, i + is + 1, n);
    recon(l, 1, 2) = reconvar(l + 0 * ndofs, k + ks, j + js - 1, i + is, n);
    recon(l, 1, 3) = reconvar(l + 0 * ndofs, k + ks, j + js - 1, i + is + 1, n);
    recon(l, 1, 4) = reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
  }
  Q2D<ndofs>(vel, recon, flux);

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l + 0 * ndofs, k + ks, j + js, i + is, n) = vel(l, 0);
      qflux(l + 1 * ndofs, k + ks, j + js, i + is, n) = vel(l, 1);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l + 0 * ndofs, k + ks, j + js, i + is, n) += vel(l, 0);
      qflux(l + 1 * ndofs, k + ks, j + js, i + is, n) += vel(l, 1);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Q_nonEC(const real5d &qflux, const real5d &reconvar,
                                 const real5d &Uvar, int is, int js, int ks,
                                 int i, int j, int k, int n) {

  SArray<real, 2, ndofs, 2> vel;
  SArray<real, 2, 2, 4> flux;
  flux(0, 0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(0, 1) = Uvar(1, k + ks, j + js, i + is - 1, n);
  flux(0, 2) = Uvar(1, k + ks, j + js + 1, i + is, n);
  flux(0, 3) = Uvar(1, k + ks, j + js + 1, i + is - 1, n);

  flux(1, 0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1, 1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  flux(1, 2) = Uvar(0, k + ks, j + js - 1, i + is, n);
  flux(1, 3) = Uvar(0, k + ks, j + js - 1, i + is + 1, n);

  SArray<real, 2, ndofs, 2> recon;
  for (int l = 0; l < ndofs; l++) {
    recon(l, 0) = reconvar(l + 0 * ndofs, k + ks, j + js, i + is, n);
    recon(l, 1) = reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
  }
  Q2D_nonEC<ndofs>(vel, recon, flux);

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l + 0 * ndofs, k + ks, j + js, i + is, n) = vel(l, 0);
      qflux(l + 1 * ndofs, k + ks, j + js, i + is, n) = vel(l, 1);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l + 0 * ndofs, k + ks, j + js, i + is, n) += vel(l, 0);
      qflux(l + 1 * ndofs, k + ks, j + js, i + is, n) += vel(l, 1);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_w_EC(const real5d &qflux, const real5d &reconvar,
                                  const real5d &vertreconvar,
                                  const real5d &Uvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 4> recon;
  flux(0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  flux(2) = Uvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = Uvar(0, k + ks + 1, j + js, i + is + 1, n);
  for (int l = 0; l < ndofs; l++) {
    recon(0) = (vertreconvar(l, k + ks, j + js, i + is, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (vertreconvar(l, k + ks, j + js, i + is + 1, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(2) = (vertreconvar(l, k + ks + 1, j + js, i + is, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(3) = (vertreconvar(l, k + ks + 1, j + js, i + is + 1, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;

    const int sgn = ndims > 1 ? -1 : 1;
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp *
          (flux(0) * recon(0) + flux(1) * recon(1) + flux(2) * recon(2) +
           flux(3) * recon(3));
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp *
          (flux(0) * recon(0) + flux(1) * recon(1) + flux(2) * recon(2) +
           flux(3) * recon(3));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_w_EC_top(const real5d &qflux,
                                      const real5d &reconvar,
                                      const real5d &vertreconvar,
                                      const real5d &Uvar, int is, int js,
                                      int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 2> recon;
  flux(0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  for (int l = 0; l < ndofs; l++) {
    recon(0) = (vertreconvar(l, k + ks, j + js, i + is, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (vertreconvar(l, k + ks, j + js, i + is + 1, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    const int sgn = ndims > 1 ? -1 : 1;
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_w_EC_bottom(const real5d &qflux,
                                         const real5d &reconvar,
                                         const real5d &vertreconvar,
                                         const real5d &Uvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 2> recon;
  flux(0) = Uvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks + 1, j + js, i + is + 1, n);
  for (int l = 0; l < ndofs; l++) {
    recon(0) = (vertreconvar(l, k + ks + 1, j + js, i + is, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (vertreconvar(l, k + ks + 1, j + js, i + is + 1, n) +
                reconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    const int sgn = ndims > 1 ? -1 : 1;
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_w_nonEC(const real5d &qflux,
                                     const real5d &reconvar, const real5d &Uvar,
                                     int is, int js, int ks, int i, int j,
                                     int k, int n) {
  SArray<real, 1, 4> flux;
  flux(0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  flux(2) = Uvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = Uvar(0, k + ks + 1, j + js, i + is + 1, n);

  const int sgn = ndims > 1 ? -1 : 1;

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          reconvar(l, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          reconvar(l, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_w_nonEC_top(const real5d &qflux,
                                         const real5d &reconvar,
                                         const real5d &Uvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  flux(0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks, j + js, i + is + 1, n);

  const int sgn = ndims > 1 ? -1 : 1;

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp * (flux(0) + flux(1)) *
          reconvar(l, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp * (flux(0) + flux(1)) *
          reconvar(l, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_w_nonEC_bottom(const real5d &qflux,
                                            const real5d &reconvar,
                                            const real5d &Uvar, int is, int js,
                                            int ks, int i, int j, int k,
                                            int n) {
  SArray<real, 1, 2> flux;
  flux(0) = Uvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks + 1, j + js, i + is + 1, n);

  const int sgn = ndims > 1 ? -1 : 1;

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp * (flux(0) + flux(1)) *
          reconvar(l, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp * (flux(0) + flux(1)) *
          reconvar(l, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_w_EC(const real5d &qflux, const real5d &reconvar,
                                  const real5d &vertreconvar,
                                  const real5d &Uvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 4> recon;
  flux(0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks, j + js + 1, i + is, n);
  flux(2) = Uvar(1, k + ks + 1, j + js, i + is, n);
  flux(3) = Uvar(1, k + ks + 1, j + js + 1, i + is, n);

  for (int l = 0; l < ndofs; l++) {
    recon(0) = (vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n) +
                reconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (vertreconvar(l + ndofs * 1, k + ks, j + js + 1, i + is, n) +
                reconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(2) = (vertreconvar(l + ndofs * 1, k + ks + 1, j + js, i + is, n) +
                reconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(3) = (vertreconvar(l + ndofs * 1, k + ks + 1, j + js + 1, i + is, n) +
                reconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;

    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1) +
                     flux(2) * recon(2) + flux(3) * recon(3));
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1) +
                     flux(2) * recon(2) + flux(3) * recon(3));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_w_EC_top(const real5d &qflux,
                                      const real5d &reconvar,
                                      const real5d &vertreconvar,
                                      const real5d &Uvar, int is, int js,
                                      int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 2> recon;
  flux(0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks, j + js + 1, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    recon(0) = (vertreconvar(l + 1 * ndofs, k + ks, j + js, i + is, n) +
                reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (vertreconvar(l + 1 * ndofs, k + ks, j + js + 1, i + is, n) +
                reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n)) *
               0.5_fp;
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_w_EC_bottom(const real5d &qflux,
                                         const real5d &reconvar,
                                         const real5d &vertreconvar,
                                         const real5d &Uvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 2> recon;
  flux(0) = Uvar(1, k + ks + 1, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks + 1, j + js + 1, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    recon(0) = (vertreconvar(l + 1 * ndofs, k + ks + 1, j + js, i + is, n) +
                reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (vertreconvar(l + 1 * ndofs, k + ks + 1, j + js + 1, i + is, n) +
                reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n)) *
               0.5_fp;
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_w_nonEC(const real5d &qflux,
                                     const real5d &reconvar, const real5d &Uvar,
                                     int is, int js, int ks, int i, int j,
                                     int k, int n) {
  SArray<real, 1, 4> flux;
  flux(0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks, j + js + 1, i + is, n);
  flux(2) = Uvar(1, k + ks + 1, j + js, i + is, n);
  flux(3) = Uvar(1, k + ks + 1, j + js + 1, i + is, n);

  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          reconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          reconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_w_nonEC_top(const real5d &qflux,
                                         const real5d &reconvar,
                                         const real5d &Uvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  flux(0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks, j + js + 1, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          0.25_fp * (flux(0) * flux(1)) *
          reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          0.25_fp * (flux(0) * flux(1)) *
          reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_w_nonEC_bottom(const real5d &qflux,
                                            const real5d &reconvar,
                                            const real5d &Uvar, int is, int js,
                                            int ks, int i, int j, int k,
                                            int n) {
  SArray<real, 1, 2> flux;
  flux(0) = Uvar(1, k + ks + 1, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks + 1, j + js + 1, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      qflux(l, k + ks, j + js, i + is, n) =
          0.25_fp * (flux(0) + flux(1)) *
          reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    }
    if (addmode == ADD_MODE::ADD) {
      qflux(l, k + ks, j + js, i + is, n) +=
          0.25_fp * (flux(0) + flux(1)) *
          reconvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_EC(const real5d &qvertflux,
                                  const real5d &reconvar,
                                  const real5d &vertreconvar,
                                  const real5d &UWvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 4> recon;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js, i + is - 1, n);
  flux(2) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = UWvar(0, k + ks + 1, j + js, i + is - 1, n);
  for (int l = 0; l < ndofs; l++) {
    // Have to subtract 1 in k here because UW has an extra dof at boundary
    // compared to v!
    recon(0) = (reconvar(l, k + ks - 1, j + js, i + is, n) +
                vertreconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (reconvar(l, k + ks - 1, j + js, i + is - 1, n) +
                vertreconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(2) = (reconvar(l, k + ks, j + js, i + is, n) +
                vertreconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(3) = (reconvar(l, k + ks, j + js, i + is - 1, n) +
                vertreconvar(l, k + ks, j + js, i + is, n)) *
               0.5_fp;

    // Added the minus sign here
    const int sgn = ndims > 1 ? 1 : -1;
    if (addmode == ADD_MODE::REPLACE) {
      qvertflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp *
          (flux(0) * recon(0) + flux(1) * recon(1) + flux(2) * recon(2) +
           flux(3) * recon(3));
    }
    if (addmode == ADD_MODE::ADD) {
      qvertflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp *
          (flux(0) * recon(0) + flux(1) * recon(1) + flux(2) * recon(2) +
           flux(3) * recon(3));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_top(const real5d &qvertflux,
                                   const real5d &vertreconvar,
                                   const real5d &UWvar, int is, int js, int ks,
                                   int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  flux(0) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks + 1, j + js, i + is - 1, n);
  // Added the minus sign here
  const int sgn = ndims > 1 ? 1 : -1;
  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l, k + ks, j + js, i + is, n);
    }
  }
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_bottom(const real5d &qvertflux,
                                      const real5d &vertreconvar,
                                      const real5d &UWvar, int is, int js,
                                      int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js, i + is - 1, n);
  // Added the minus sign here
  const int sgn = ndims > 1 ? 1 : -1;
  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l, k + ks, j + js, i + is, n);
    }
  }
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_nonEC(const real5d &qvertflux,
                                     const real5d &vertreconvar,
                                     const real5d &UWvar, int is, int js,
                                     int ks, int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js, i + is - 1, n);
  flux(2) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = UWvar(0, k + ks + 1, j + js, i + is - 1, n);
  for (int l = 0; l < ndofs; l++) {
    // Added the minus sign here
    const int sgn = ndims > 1 ? 1 : -1;
    if (addmode == ADD_MODE::REPLACE) {
      qvertflux(l, k + ks, j + js, i + is, n) =
          sgn * 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          vertreconvar(l, k + ks, j + js, i + is, n);
    }
    if (addmode == ADD_MODE::ADD) {
      qvertflux(l, k + ks, j + js, i + is, n) +=
          sgn * 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          vertreconvar(l, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_nonEC_top(const real5d &qvertflux,
                                         const real5d &vertreconvar,
                                         const real5d &UWvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  compute_Qxz_u_top<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks,
                                    i, j, k, n);
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_EC_top(const real5d &qvertflux,
                                      const real5d &reconvar,
                                      const real5d &vertreconvar,
                                      const real5d &UWvar, int is, int js,
                                      int ks, int i, int j, int k, int n) {
  compute_Qxz_u_top<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks,
                                    i, j, k, n);
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_nonEC_bottom(const real5d &qvertflux,
                                            const real5d &vertreconvar,
                                            const real5d &UWvar, int is, int js,
                                            int ks, int i, int j, int k,
                                            int n) {
  compute_Qxz_u_bottom<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js,
                                       ks, i, j, k, n);
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qxz_u_EC_bottom(const real5d &qvertflux,
                                         const real5d &reconvar,
                                         const real5d &vertreconvar,
                                         const real5d &UWvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  compute_Qxz_u_bottom<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js,
                                       ks, i, j, k, n);
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_EC(const real5d &qvertflux,
                                  const real5d &reconvar,
                                  const real5d &vertreconvar,
                                  const real5d &UWvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 4> recon;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js - 1, i + is, n);
  flux(2) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = UWvar(0, k + ks + 1, j + js - 1, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    // Have to subtract 1 in k here because UW has an extra dof at boundary
    // compared to v!
    recon(0) = (reconvar(l + ndofs * 1, k + ks - 1, j + js, i + is, n) +
                vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(1) = (reconvar(l + ndofs * 1, k + ks - 1, j + js - 1, i + is, n) +
                vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(2) = (reconvar(l + ndofs * 1, k + ks, j + js, i + is, n) +
                vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    recon(3) = (reconvar(l + ndofs * 1, k + ks, j + js - 1, i + is, n) +
                vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n)) *
               0.5_fp;
    if (addmode == ADD_MODE::REPLACE) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) =
          -0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1) +
                      flux(2) * recon(2) + flux(3) * recon(3));
    }
    if (addmode == ADD_MODE::ADD) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) +=
          -0.25_fp * (flux(0) * recon(0) + flux(1) * recon(1) +
                      flux(2) * recon(2) + flux(3) * recon(3));
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_top(const real5d &qvertflux,
                                   const real5d &vertreconvar,
                                   const real5d &UWvar, int is, int js, int ks,
                                   int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  flux(0) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks + 1, j + js - 1, i + is, n);
  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) =
          -0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) +=
          -0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
  }
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_bottom(const real5d &qvertflux,
                                      const real5d &vertreconvar,
                                      const real5d &UWvar, int is, int js,
                                      int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> flux;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js - 1, i + is, n);
  // Added the minus sign here
  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) =
          -0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) +=
          -0.5_fp * (flux(0) + flux(1)) *
          vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
  }
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_nonEC(const real5d &qvertflux,
                                     const real5d &vertreconvar,
                                     const real5d &UWvar, int is, int js,
                                     int ks, int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js - 1, i + is, n);
  flux(2) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = UWvar(0, k + ks + 1, j + js - 1, i + is, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) =
          -0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
    if (addmode == ADD_MODE::ADD) {
      qvertflux(l + ndofs * 1, k + ks, j + js, i + is, n) +=
          -0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3)) *
          vertreconvar(l + ndofs * 1, k + ks, j + js, i + is, n);
    }
  }
}

template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_nonEC_top(const real5d &qvertflux,
                                         const real5d &vertreconvar,
                                         const real5d &UWvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  compute_Qyz_v_top<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks,
                                    i, j, k, n);
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_EC_top(const real5d &qvertflux,
                                      const real5d &reconvar,
                                      const real5d &vertreconvar,
                                      const real5d &UWvar, int is, int js,
                                      int ks, int i, int j, int k, int n) {
  compute_Qyz_v_top<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js, ks,
                                    i, j, k, n);
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_nonEC_bottom(const real5d &qvertflux,
                                            const real5d &vertreconvar,
                                            const real5d &UWvar, int is, int js,
                                            int ks, int i, int j, int k,
                                            int n) {
  compute_Qyz_v_bottom<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js,
                                       ks, i, j, k, n);
}
template <index_t ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_Qyz_v_EC_bottom(const real5d &qvertflux,
                                         const real5d &reconvar,
                                         const real5d &vertreconvar,
                                         const real5d &UWvar, int is, int js,
                                         int ks, int i, int j, int k, int n) {
  compute_Qyz_v_bottom<ndofs, addmode>(qvertflux, vertreconvar, UWvar, is, js,
                                       ks, i, j, k, n);
}

// W
void YAKL_INLINE W2D(SArray<real, 1, 2> &vel,
                     SArray<real, 2, 2, 4> const &flux) {
  // x-dir
  vel(0) = -0.25_fp * (flux(0, 0) + flux(0, 1) + flux(0, 2) + flux(0, 3));
  // y-dir
  vel(1) = 0.25_fp * (flux(1, 0) + flux(1, 1) + flux(1, 2) + flux(1, 3));
}

void YAKL_INLINE compute_W(const real5d &UTvar, const real5d &Uvar, int is,
                           int js, int ks, int i, int j, int k, int n) {

  SArray<real, 1, 2> ut;
  SArray<real, 2, 2, 4> flux;
  flux(0, 0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(0, 1) = Uvar(1, k + ks, j + js, i + is - 1, n);
  flux(0, 2) = Uvar(1, k + ks, j + js + 1, i + is, n);
  flux(0, 3) = Uvar(1, k + ks, j + js + 1, i + is - 1, n);

  flux(1, 0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1, 1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  flux(1, 2) = Uvar(0, k + ks, j + js - 1, i + is, n);
  flux(1, 3) = Uvar(0, k + ks, j + js - 1, i + is + 1, n);

  W2D(ut, flux);
  UTvar(0, k + ks, j + js, i + is, n) = ut(0);
  UTvar(1, k + ks, j + js, i + is, n) = ut(1);
}

// Wxz
void YAKL_INLINE Wxz_u(SArray<real, 1, 1> &vel,
                       SArray<real, 1, 4> const &flux) {
  // Added the minus sign here
  vel(0) = -0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3));
}
void YAKL_INLINE Wxz_u_boundary(SArray<real, 1, 1> &vel,
                                SArray<real, 1, 2> const &flux) {

  // Added the minus sign here
  vel(0) = -0.5_fp * (flux(0) + flux(1));
}

void YAKL_INLINE compute_Wxz_u(const real5d &VTvar, const real5d &UWvar, int is,
                               int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 1> vel;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js, i + is - 1, n);
  flux(2) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = UWvar(0, k + ks + 1, j + js, i + is - 1, n);

  Wxz_u(vel, flux);
  VTvar(0, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wxz_u_top(const real5d &VTvar, const real5d &UWvar,
                                   int is, int js, int ks, int i, int j, int k,
                                   int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks + 1, j + js, i + is - 1, n);

  Wxz_u_boundary(vel, flux);
  VTvar(0, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wxz_u_bottom(const real5d &VTvar, const real5d &UWvar,
                                      int is, int js, int ks, int i, int j,
                                      int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js, i + is - 1, n);

  Wxz_u_boundary(vel, flux);
  VTvar(0, k + ks, j + js, i + is, n) = vel(0);
}

void YAKL_INLINE Wxz_w(SArray<real, 1, 1> &vel,
                       SArray<real, 1, 4> const &flux) {
  vel(0) = 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3));
}
void YAKL_INLINE Wxz_w_boundary(SArray<real, 1, 1> &vel,
                                SArray<real, 1, 2> const &flux) {

  vel(0) = 0.25_fp * (flux(0) + flux(1));
}

void YAKL_INLINE compute_Wxz_w(const real5d &WTvar, const real5d &Uvar, int is,
                               int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 1> vel;
  flux(0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks, j + js, i + is + 1, n);
  flux(2) = Uvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = Uvar(0, k + ks + 1, j + js, i + is + 1, n);

  Wxz_w(vel, flux);
  WTvar(0, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wxz_w_top(const real5d &WTvar, const real5d &Uvar,
                                   int is, int js, int ks, int i, int j, int k,
                                   int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = Uvar(0, k + ks, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks, j + js, i + is + 1, n);

  Wxz_w_boundary(vel, flux);
  WTvar(0, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wxz_w_bottom(const real5d &WTvar, const real5d &Uvar,
                                      int is, int js, int ks, int i, int j,
                                      int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = Uvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = Uvar(0, k + ks + 1, j + js, i + is + 1, n);

  Wxz_w_boundary(vel, flux);
  WTvar(0, k + ks, j + js, i + is, n) = vel(0);
}

// Wyz
void YAKL_INLINE Wyz_u(SArray<real, 1, 1> &vel,
                       SArray<real, 1, 4> const &flux) {
  vel(0) = 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3));
}
void YAKL_INLINE Wyz_u_boundary(SArray<real, 1, 1> &vel,
                                SArray<real, 1, 2> const &flux) {

  vel(0) = 0.5_fp * (flux(0) + flux(1));
}
void YAKL_INLINE compute_Wyz_u(const real5d &VTvar, const real5d &UWvar, int is,
                               int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 1> vel;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js - 1, i + is, n);
  flux(2) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(3) = UWvar(0, k + ks + 1, j + js - 1, i + is, n);

  Wyz_u(vel, flux);
  VTvar(1, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wyz_u_top(const real5d &VTvar, const real5d &UWvar,
                                   int is, int js, int ks, int i, int j, int k,
                                   int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = UWvar(0, k + ks + 1, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks + 1, j + js - 1, i + is, n);

  Wyz_u_boundary(vel, flux);
  VTvar(1, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wyz_u_bottom(const real5d &VTvar, const real5d &UWvar,
                                      int is, int js, int ks, int i, int j,
                                      int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = UWvar(0, k + ks, j + js, i + is, n);
  flux(1) = UWvar(0, k + ks, j + js - 1, i + is, n);

  Wyz_u_boundary(vel, flux);
  VTvar(1, k + ks, j + js, i + is, n) = vel(0);
}

void YAKL_INLINE Wyz_w(SArray<real, 1, 1> &vel,
                       SArray<real, 1, 4> const &flux) {
  vel(0) = 0.25_fp * (flux(0) + flux(1) + flux(2) + flux(3));
}
void YAKL_INLINE Wyz_w_boundary(SArray<real, 1, 1> &vel,
                                SArray<real, 1, 2> const &flux) {

  vel(0) = 0.25_fp * (flux(0) + flux(1));
}

void YAKL_INLINE compute_Wyz_w(const real5d &WTvar, const real5d &Uvar, int is,
                               int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, 4> flux;
  SArray<real, 1, 1> vel;
  flux(0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks, j + js + 1, i + is, n);
  flux(2) = Uvar(1, k + ks + 1, j + js, i + is, n);
  flux(3) = Uvar(1, k + ks + 1, j + js + 1, i + is, n);

  Wyz_w(vel, flux);
  WTvar(1, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wyz_w_top(const real5d &WTvar, const real5d &Uvar,
                                   int is, int js, int ks, int i, int j, int k,
                                   int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = Uvar(1, k + ks, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks, j + js + 1, i + is, n);

  Wyz_w_boundary(vel, flux);
  WTvar(1, k + ks, j + js, i + is, n) = vel(0);
}
void YAKL_INLINE compute_Wyz_w_bottom(const real5d &WTvar, const real5d &Uvar,
                                      int is, int js, int ks, int i, int j,
                                      int k, int n) {
  SArray<real, 1, 2> flux;
  SArray<real, 1, 1> vel;
  flux(0) = Uvar(1, k + ks + 1, j + js, i + is, n);
  flux(1) = Uvar(1, k + ks + 1, j + js + 1, i + is, n);

  Wyz_w_boundary(vel, flux);
  WTvar(1, k + ks, j + js, i + is, n) = vel(0);
}

// R
void YAKL_INLINE R(SArray<real, 1, 1> &vard, SArray<real, 1, 2> const &varp) {
  vard(0) = 0.5_fp * (varp(0) + varp(1));
}
void YAKL_INLINE R(SArray<real, 1, 1> &vard, SArray<real, 1, 4> const &varp) {
  vard(0) = 0.25_fp * (varp(0) + varp(1) + varp(2) + varp(3));
}
void YAKL_INLINE R(SArray<real, 1, 1> &vard, SArray<real, 1, 8> const &varp) {
  vard(0) = 0.125_fp * (varp(0) + varp(1) + varp(2) + varp(3) + varp(4) +
                        varp(5) + varp(6) + varp(7));
}

template <index_t dof>
void YAKL_INLINE compute_R(real &var, const real5d &densvar, int is, int js,
                           int ks, int i, int j, int k, int n) {
  SArray<real, 1, 1> dualdens;
  if (ndims == 1) {
    SArray<real, 1, 2> dens;
    dens(0) = densvar(dof, k + ks, j + js, i + is, n);
    dens(1) = densvar(dof, k + ks, j + js, i + is - 1, n);
    R(dualdens, dens);
  }
  if (ndims == 2) {
    SArray<real, 1, 4> dens;
    dens(0) = densvar(dof, k + ks, j + js, i + is, n);
    dens(1) = densvar(dof, k + ks, j + js, i + is - 1, n);
    dens(2) = densvar(dof, k + ks, j + js - 1, i + is, n);
    dens(3) = densvar(dof, k + ks, j + js - 1, i + is - 1, n);
    R(dualdens, dens);
  }
  if (ndims == 3) {
    SArray<real, 1, 8> dens;
    dens(0) = densvar(dof, k + ks, j + js, i + is, n);
    dens(1) = densvar(dof, k + ks, j + js, i + is - 1, n);
    dens(2) = densvar(dof, k + ks, j + js - 1, i + is, n);
    dens(3) = densvar(dof, k + ks, j + js - 1, i + is - 1, n);
    dens(4) = densvar(dof, k + ks - 1, j + js, i + is, n);
    dens(5) = densvar(dof, k + ks - 1, j + js, i + is - 1, n);
    dens(6) = densvar(dof, k + ks - 1, j + js - 1, i + is, n);
    dens(7) = densvar(dof, k + ks - 1, j + js - 1, i + is - 1, n);
    R(dualdens, dens);
  }
  var = dualdens(0);
}

template <index_t dof>
void YAKL_INLINE compute_R(const real5d &var, const real5d &densvar, int is,
                           int js, int ks, int i, int j, int k, int n) {
  real var2;
  compute_R<dof>(var2, densvar, is, js, ks, i, j, k, n);
  var(dof, k + ks, j + js, i + is, n) = var2;
}

void YAKL_INLINE Rbnd(SArray<real, 1, 1> &vard,
                      SArray<real, 1, 4> const &varp) {
  vard(0) = 0.25_fp * (varp(0) + varp(1)) + 0.5_fp * (varp(2) + varp(3));
}

// phi
void YAKL_INLINE phi(SArray<real, 1, ndims> &vare,
                     SArray<real, 2, ndims, 2> const &var0) {
  for (int d = 0; d < ndims; d++) {
    vare(d) = 0.5_fp * (var0(d, 0) + var0(d, 1));
  }
}

void YAKL_INLINE compute_phi(const real5d &var, const real5d &densvar, int is,
                             int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndims> xe;
  SArray<real, 2, ndims, 2> x0;
  for (int d = 0; d < ndims; d++) {
    if (d == 0) {
      x0(d, 0) = densvar(0, k + ks, j + js, i + is, n);
      x0(d, 1) = densvar(0, k + ks, j + js, i + is - 1, n);
    }
    if (d == 1) {
      x0(d, 0) = densvar(0, k + ks, j + js, i + is, n);
      x0(d, 1) = densvar(0, k + ks, j + js - 1, i + is, n);
    }
    if (d == 2) {
      x0(d, 0) = densvar(0, k + ks, j + js, i + is, n);
      x0(d, 1) = densvar(0, k + ks - 1, j + js, i + is, n);
    }
  }
  phi(xe, x0);
  for (int d = 0; d < ndims; d++) {
    var(d, k + ks, j + js, i + is, n) = xe(d);
  }
}

real YAKL_INLINE phiW(SArray<real, 1, 2> const &var0) {
  return 0.5_fp * (var0(0) + var0(1));
}

void YAKL_INLINE compute_phiW(const real5d &var, const real5d &densvar, int is,
                              int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, 2> x0;
  x0(0) = densvar(0, k + ks, j + js, i + is, n);
  x0(1) = densvar(0, k + ks - 1, j + js, i + is, n);
  var(0, k + ks, j + js, i + is, n) = phiW(x0);
}

// phiT
void YAKL_INLINE phiT(SArray<real, 1, 1> &ke,
                      SArray<real, 2, ndims, 2> const &u,
                      SArray<real, 2, ndims, 2> const &v) {
  ke(0) = 0.0;
  for (int d = 0; d < ndims; d++) {
    ke(0) += v(d, 0) * u(d, 0) + v(d, 1) * u(d, 1);
  }
  ke(0) *= 0.5_fp;
}

template <ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_phiT(const real5d &var, const real5d &uvar,
                              const real5d &vvar, int is, int js, int ks, int i,
                              int j, int k, int n) {
  SArray<real, 1, 1> ke;
  SArray<real, 2, ndims, 2> u;
  SArray<real, 2, ndims, 2> v;
  for (int d = 0; d < ndims; d++) {
    if (d == 0) {
      u(d, 0) = uvar(d, k + ks, j + js, i + is, n);
      u(d, 1) = uvar(d, k + ks, j + js, i + is + 1, n);
      v(d, 0) = vvar(d, k + ks, j + js, i + is, n);
      v(d, 1) = vvar(d, k + ks, j + js, i + is + 1, n);
    }
    if (d == 1) {
      u(d, 0) = uvar(d, k + ks, j + js, i + is, n);
      u(d, 1) = uvar(d, k + ks, j + js + 1, i + is, n);
      v(d, 0) = vvar(d, k + ks, j + js, i + is, n);
      v(d, 1) = vvar(d, k + ks, j + js + 1, i + is, n);
    }
    if (d == 2) {
      u(d, 0) = uvar(d, k + ks, j + js, i + is, n);
      u(d, 1) = uvar(d, k + ks + 1, j + js, i + is, n);
      v(d, 0) = vvar(d, k + ks, j + js, i + is, n);
      v(d, 1) = vvar(d, k + ks + 1, j + js, i + is, n);
    }
  }
  phiT(ke, u, v);
  if (addmode == ADD_MODE::REPLACE) {
    var(0, k + ks, j + js, i + is, n) = ke(0);
  }
  if (addmode == ADD_MODE::ADD) {
    var(0, k + ks, j + js, i + is, n) += ke(0);
  }
}
void YAKL_INLINE compute_phiT(SArray<real, 1, 1> &var,
                              SArray<real, 2, ndims, 2> const &u,
                              const real5d &vvar, int is, int js, int ks, int i,
                              int j, int k, int n) {
  SArray<real, 2, ndims, 2> v;
  for (int d = 0; d < ndims; d++) {
    if (d == 0) {
      v(d, 0) = vvar(d, k + ks, j + js, i + is, n);
      v(d, 1) = vvar(d, k + ks, j + js, i + is + 1, n);
    }
    if (d == 1) {
      v(d, 0) = vvar(d, k + ks, j + js, i + is, n);
      v(d, 1) = vvar(d, k + ks, j + js + 1, i + is, n);
    }
    if (d == 2) {
      v(d, 0) = vvar(d, k + ks, j + js, i + is, n);
      v(d, 1) = vvar(d, k + ks + 1, j + js, i + is, n);
    }
  }
  phiT(var, u, v);
}

real YAKL_INLINE phiTW(SArray<real, 1, 2> const &u,
                       SArray<real, 1, 2> const &v) {
  return (v(0) * u(0) + v(1) * u(1)) * 0.5_fp;
}

template <ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_phiTW(const real5d &var, const real5d &uwvar,
                               const real5d &wvar, int is, int js, int ks,
                               int i, int j, int k, int n) {
  SArray<real, 1, 2> u;
  SArray<real, 1, 2> v;
  u(0) = uwvar(0, k + ks, j + js, i + is, n);
  u(1) = uwvar(0, k + ks + 1, j + js, i + is, n);
  // Have to subtract 1 from k here since UW has an extra dof compared to w
  v(0) = wvar(0, k + ks - 1, j + js, i + is, n);
  v(1) = wvar(0, k + ks, j + js, i + is, n);
  if (addmode == ADD_MODE::REPLACE) {
    var(0, k + ks, j + js, i + is, n) = phiTW(u, v);
  }
  if (addmode == ADD_MODE::ADD) {
    var(0, k + ks, j + js, i + is, n) += phiTW(u, v);
  }
}

template <ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_phiTW(SArray<real, 1, 1> &var,
                               SArray<real, 1, 2> const &uw, const real5d &wvar,
                               int is, int js, int ks, int i, int j, int k,
                               int n) {
  SArray<real, 1, 2> v;
  // SHOULD SUBTRACT HERE ALSO?
  v(0) = wvar(0, k + ks - 1, j + js, i + is, n);
  v(1) = wvar(0, k + ks, j + js, i + is, n);
  if (addmode == ADD_MODE::REPLACE) {
    var(0) = phiTW(uw, v);
  }
  if (addmode == ADD_MODE::ADD) {
    var(0) += phiTW(uw, v);
  }
}
} // namespace pamc
