#pragma once

#include "common.h"

template <uint ndofs>
void YAKL_INLINE D0(SArray<real, 1, ndims> &var,
                    SArray<real, 3, ndofs, ndims, 2> const &dens) {

  for (int d = 0; d < ndims; d++) {
    var(d) = 0.0;
    for (int l = 0; l < ndofs; l++) {
      var(d) += (dens(l, d, 1) - dens(l, d, 0));
    }
  }
}

template <uint ndofs>
void YAKL_INLINE compute_D0(SArray<real, 1, ndims> &tend, const real5d &densvar,
                            int is, int js, int ks, int i, int j, int k,
                            int n) {
  SArray<real, 3, ndofs, ndims, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < ndims; d++) {
      dens(l, d, 1) = densvar(l, k + ks, j + js, i + is, n);
      if (d == 0) {
        dens(l, d, 0) = densvar(l, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        dens(l, d, 0) = densvar(l, k + ks, j + js - 1, i + is, n);
      }
      // if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }
  }
  D0<ndofs>(tend, dens);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D0(const real5d &tendvar, const real5d &densvar,
                            int is, int js, int ks, int i, int j, int k,
                            int n) {
  SArray<real, 1, ndims> tend;
  compute_D0<ndofs>(tend, densvar, is, js, ks, i, j, k, n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(d, k + ks, j + js, i + is, n) = tend(d);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(d, k + ks, j + js, i + is, n) += tend(d);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE wD0(SArray<real, 1, ndims> &var,
                     SArray<real, 2, ndofs, ndims> const &recon,
                     SArray<real, 3, ndofs, ndims, 2> const &dens) {

  for (int d = 0; d < ndims; d++) {
    var(d) = 0.0;
    for (int l = 0; l < ndofs; l++) {
      var(d) += recon(l, d) * (dens(l, d, 1) - dens(l, d, 0));
    }
  }
}

template <uint ndofs, class R>
void YAKL_INLINE compute_wD0(SArray<real, 1, ndims> &tend, const R &reconvar,
                             const real5d &densvar, int is, int js, int ks,
                             int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < ndims; d++) {
      dens(l, d, 1) = densvar(l, k + ks, j + js, i + is, n);
      if (d == 0) {
        dens(l, d, 0) = densvar(l, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        dens(l, d, 0) = densvar(l, k + ks, j + js - 1, i + is, n);
      }
      // if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }
  }

  // we need to load recon
  if constexpr (std::is_same_v<R, real5d> || std::is_same_v<R, real3d>) {
    SArray<real, 2, ndofs, ndims> recon;
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        // full state recon
        if constexpr (std::is_same_v<R, real5d>) {
          recon(l, d) = reconvar(d + l * ndims, k + ks, j + js, i + is, n);
        }
        // reference state
        if constexpr (std::is_same_v<R, real3d>) {
          recon(l, d) = reconvar(d + l * ndims, k + ks, n);
        }
      }
    }
    wD0<ndofs>(tend, recon, dens);
  } else {
    wD0<ndofs>(tend, reconvar, dens);
  }
}

template <uint ndofs, class R>
void YAKL_INLINE compute_wD0(SArray<real, 1, ndims> &tend, const R &reconvar,
                             const SArray<int, 1, ndofs> &active_dens_ids,
                             const real5d &densvar, int is, int js, int ks,
                             int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < ndims; d++) {
      dens(l, d, 1) = densvar(l, k + ks, j + js, i + is, n);
      if (d == 0) {
        dens(l, d, 0) = densvar(l, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        dens(l, d, 0) = densvar(l, k + ks, j + js - 1, i + is, n);
      }
      // if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }
  }

  // we need to load recon
  if constexpr (std::is_same_v<R, real5d> || std::is_same_v<R, real3d>) {
    SArray<real, 2, ndofs, ndims> recon;
    for (int l = 0; l < ndofs; l++) {
      int al = active_dens_ids(l);
      for (int d = 0; d < ndims; d++) {
        // full state recon
        if constexpr (std::is_same_v<R, real5d>) {
          recon(l, d) = reconvar(d + al * ndims, k + ks, j + js, i + is, n);
        }
        // reference state
        if constexpr (std::is_same_v<R, real3d>) {
          recon(l, d) = reconvar(d + al * ndims, k + ks, n);
        }
      }
    }
    wD0<ndofs>(tend, recon, dens);
  } else {
    wD0<ndofs>(tend, reconvar, dens);
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
void YAKL_INLINE compute_wD0(const real5d &tendvar, const R &reconvar,
                             const real5d &densvar, int is, int js, int ks,
                             int i, int j, int k, int n) {
  SArray<real, 1, ndims> tend;
  compute_wD0<ndofs>(tend, reconvar, densvar, is, js, ks, i, j, k, n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(d, k + ks, j + js, i + is, n) = tend(d);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(d, k + ks, j + js, i + is, n) += tend(d);
    }
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
void YAKL_INLINE compute_wD0(const real5d &tendvar, const R &reconvar,
                             const SArray<int, 1, ndofs> &active_dens_ids,
                             const real5d &densvar, int is, int js, int ks,
                             int i, int j, int k, int n) {
  SArray<real, 1, ndims> tend;
  compute_wD0<ndofs>(tend, reconvar, active_dens_ids, densvar, is, js, ks, i, j,
                     k, n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(d, k + ks, j + js, i + is, n) = tend(d);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(d, k + ks, j + js, i + is, n) += tend(d);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE D0_vert(real &var, SArray<real, 2, ndofs, 2> const &dens) {

  var = 0.0;
  for (int l = 0; l < ndofs; l++) {
    var += (dens(l, 1) - dens(l, 0));
  }
}

template <uint ndofs>
real YAKL_INLINE compute_D0_vert(const real5d &densvar, int is, int js, int ks,
                                 int i, int j, int k, int n) {
  real tend;
  SArray<real, 2, ndofs, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    dens(l, 1) = densvar(l, k + ks + 1, j + js, i + is, n);
    dens(l, 0) = densvar(l, k + ks, j + js, i + is, n);
  }

  D0_vert<ndofs>(tend, dens);
  return tend;
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D0_vert(const real5d &tendvar, const real5d &densvar,
                                 int is, int js, int ks, int i, int j, int k,
                                 int n) {
  real tend = compute_D0_vert<ndofs>(densvar, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    tendvar(0, k + ks, j + js, i + is, n) = tend;
  }
  if (addmode == ADD_MODE::ADD) {
    tendvar(0, k + ks, j + js, i + is, n) += tend;
  }
}

template <uint ndofs>
void YAKL_INLINE wD0_vert(real &var, SArray<real, 1, ndofs> const &recon,
                          SArray<real, 2, ndofs, 2> const &dens) {

  var = 0.0;
  for (int l = 0; l < ndofs; l++) {
    var += recon(l) * (dens(l, 1) - dens(l, 0));
  }
}

template <uint ndofs, class R, class D>
real YAKL_INLINE compute_wD0_vert(const R &vertreconvar, const D &densvar,
                                  int is, int js, int ks, int i, int j, int k,
                                  int n) {
  real tend;
  SArray<real, 2, ndofs, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    // full state dens
    if constexpr (std::is_same_v<D, real5d>) {
      dens(l, 1) = densvar(l, k + ks + 1, j + js, i + is, n);
      dens(l, 0) = densvar(l, k + ks, j + js, i + is, n);
    }
    // reference state dens
    if constexpr (std::is_same_v<D, real3d>) {
      dens(l, 1) = densvar(l, k + ks + 1, n);
      dens(l, 0) = densvar(l, k + ks, n);
    }
  }

  // we have to load recon
  if constexpr (std::is_same_v<R, real5d> || std::is_same_v<R, real3d>) {
    SArray<real, 1, ndofs> recon;
    for (int l = 0; l < ndofs; l++) {
      // Need to add 1 to k here because UW has an extra dof at the bottom

      // full state recon
      if constexpr (std::is_same_v<R, real5d>) {
        recon(l) = vertreconvar(l, k + ks + 1, j + js, i + is, n);
      }

      // reference state recon
      if constexpr (std::is_same_v<R, real3d>) {
        recon(l) = vertreconvar(l, k + ks + 1, n);
      }
    }
    wD0_vert<ndofs>(tend, recon, dens);
  } else {
    wD0_vert<ndofs>(tend, vertreconvar, dens);
  }
  return tend;
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R, class D>
void YAKL_INLINE compute_wD0_vert(const real5d &tendvar, const R &vertreconvar,
                                  const D &densvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  real tend =
      compute_wD0_vert<ndofs>(vertreconvar, densvar, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    tendvar(0, k + ks, j + js, i + is, n) = tend;
  }
  if (addmode == ADD_MODE::ADD) {
    tendvar(0, k + ks, j + js, i + is, n) += tend;
  }
}

template <uint ndofs, class R, class D>
real YAKL_INLINE compute_wD0_vert(const R &vertreconvar,
                                  const SArray<int, 1, ndofs> &active_dens_ids,
                                  const D &densvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  real tend;
  SArray<real, 2, ndofs, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    // full state dens
    if constexpr (std::is_same_v<D, real5d>) {
      dens(l, 1) = densvar(l, k + ks + 1, j + js, i + is, n);
      dens(l, 0) = densvar(l, k + ks, j + js, i + is, n);
    }
    // reference state dens
    if constexpr (std::is_same_v<D, real3d>) {
      dens(l, 1) = densvar(l, k + ks + 1, n);
      dens(l, 0) = densvar(l, k + ks, n);
    }
  }

  // we have to load recon
  if constexpr (std::is_same_v<R, real5d> || std::is_same_v<R, real3d>) {
    SArray<real, 1, ndofs> recon;
    for (int l = 0; l < ndofs; l++) {
      int al = active_dens_ids(l);
      // Need to add 1 to k here because UW has an extra dof at the bottom

      // full state recon
      if constexpr (std::is_same_v<R, real5d>) {
        recon(l) = vertreconvar(al, k + ks + 1, j + js, i + is, n);
      }

      // reference state recon
      if constexpr (std::is_same_v<R, real3d>) {
        recon(l) = vertreconvar(al, k + ks + 1, n);
      }
    }
    wD0_vert<ndofs>(tend, recon, dens);
  } else {
    wD0_vert<ndofs>(tend, vertreconvar, dens);
  }
  return tend;
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R, class D>
void YAKL_INLINE compute_wD0_vert(const real5d &tendvar, const R &vertreconvar,
                                  const SArray<int, 1, ndofs> &active_dens_ids,
                                  const D &densvar, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  real tend = compute_wD0_vert<ndofs, R, D>(vertreconvar, active_dens_ids,
                                            densvar, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    tendvar(0, k + ks, j + js, i + is, n) = tend;
  }
  if (addmode == ADD_MODE::ADD) {
    tendvar(0, k + ks, j + js, i + is, n) += tend;
  }
}

template <uint ndofs>
void YAKL_INLINE D0bar(SArray<real, 1, ndims> &var,
                       SArray<real, 3, ndofs, ndims, 2> const &dens) {

  for (int d = 0; d < ndims; d++) {
    var(d) = 0.0;
    for (int l = 0; l < ndofs; l++) {
      if (d == 0) {
        var(d) += (dens(l, d, 1) - dens(l, d, 0));
      }
      if (d == 1) {
        var(d) -= (dens(l, d, 1) - dens(l, d, 0));
      }
    }
  }
}

template <uint ndofs>
void YAKL_INLINE compute_D0bar(SArray<real, 1, ndims> &tend,
                               const real5d &densvar, int is, int js, int ks,
                               int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < ndims; d++) {
      dens(l, d, 0) = densvar(l, k + ks, j + js, i + is, n);
      if (d == 0) {
        dens(l, d, 1) = densvar(l, k + ks, j + js + 1, i + is, n);
      }
      if (d == 1) {
        dens(l, d, 1) = densvar(l, k + ks, j + js, i + is + 1, n);
      }
      // if (d==2) {dens(l,d,0) = densvar(l,k+ks-1,j+js,i+is);}
    }
  }
  D0bar<ndofs>(tend, dens);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D0bar(const real5d &tendvar, const real5d &densvar,
                               int is, int js, int ks, int i, int j, int k,
                               int n) {
  SArray<real, 1, ndims> tend;
  compute_D0bar<ndofs>(tend, densvar, is, js, ks, i, j, k, n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(d, k + ks, j + js, i + is, n) = tend(d);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(d, k + ks, j + js, i + is, n) += tend(d);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE compute_D0bar_ext(SArray<real, 1, ndims> &tend,
                                   const real5d &densvar, int is, int js,
                                   int ks, int i, int j, int k, int n) {
  SArray<real, 3, ndofs, ndims, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    for (int d = 0; d < ndims; d++) {
      dens(l, d, 0) = densvar(l, k + ks, j + js, i + is + 1, n);
      if (d == 0) {
        dens(l, d, 1) = densvar(l, k + ks, j + js, i + is, n);
      }
    }
  }
  D0bar<ndofs>(tend, dens);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D0bar_ext(const real5d &tendvar, const real5d &densvar,
                                   int is, int js, int ks, int i, int j, int k,
                                   int n) {
  SArray<real, 1, ndims> tend;
  compute_D0bar_ext<ndofs>(tend, densvar, is, js, ks, i, j, k, n);
  for (int d = 0; d < ndims; d++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(d, k + ks, j + js, i + is, n) = tend(d);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(d, k + ks, j + js, i + is, n) += tend(d);
    }
  }
}
template <uint ndofs>
void YAKL_INLINE D0bar_vert(real &var, SArray<real, 2, ndofs, 2> const &dens) {

  var = 0.0;
  for (int l = 0; l < ndofs; l++) {
    var += (dens(l, 1) - dens(l, 0));
  }
}

template <uint ndofs>
real YAKL_INLINE compute_D0bar_vert(const real5d &densvar, int is, int js,
                                    int ks, int i, int j, int k, int n) {
  real tend;
  SArray<real, 2, ndofs, 2> dens;
  for (int l = 0; l < ndofs; l++) {
    dens(l, 1) = densvar(l, k + ks + 1, j + js, i + is, n);
    dens(l, 0) = densvar(l, k + ks, j + js, i + is, n);
  }

  D0bar_vert<ndofs>(tend, dens);
  return tend;
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D0bar_vert(const real5d &tendvar,
                                    const real5d &densvar, int is, int js,
                                    int ks, int i, int j, int k, int n) {
  real tend = compute_D0bar_vert<ndofs>(densvar, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    tendvar(0, k + ks, j + js, i + is, n) = tend;
  }
  if (addmode == ADD_MODE::ADD) {
    tendvar(0, k + ks, j + js, i + is, n) += tend;
  }
}

template <uint ndofs>
void YAKL_INLINE Dnm1bar(SArray<real, 1, ndofs> &var,
                         SArray<real, 2, ndims, 2> const &flux) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = 0.;
    for (int d = 0; d < ndims; d++) {
      var(l) += flux(d, 1) - flux(d, 0);
    }
  }
}

template <uint ndofs>
YAKL_INLINE void compute_Dnm1bar(SArray<real, 1, ndofs> &tend, const real5d &U,
                                 int is, int js, int ks, int i, int j, int k,
                                 int n) {
  SArray<real, 2, ndims, 2> flux;

  for (int d = 0; d < ndims; d++) {
    for (int m = 0; m < 2; m++) {
      if (d == 0) {
        flux(d, m) = U(d, k + ks, j + js, i + is + m, n);
      }
      if (d == 1) {
        flux(d, m) = U(d, k + ks, j + js + m, i + is, n);
      }
    }
  }
  Dnm1bar<ndofs>(tend, flux);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
YAKL_INLINE void compute_Dnm1bar(const real5d &tendvar, const real5d &U, int is,
                                 int js, int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_Dnm1bar<ndofs>(tend, U, is, js, ks, i, j, k, n);

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE wDnm1bar(SArray<real, 1, ndofs> &var,
                          SArray<real, 3, ndofs, ndims, 2> const &recon,
                          SArray<real, 2, ndims, 2> const &flux) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = 0.;
    for (int d = 0; d < ndims; d++) {
      var(l) += flux(d, 1) * recon(l, d, 1) - flux(d, 0) * recon(l, d, 0);
    }
  }
}

template <uint ndofs, class R>
YAKL_INLINE void compute_wDnm1bar(SArray<real, 1, ndofs> &tend,
                                  const R &reconvar, const real5d &U, int is,
                                  int js, int ks, int i, int j, int k, int n) {
  SArray<real, 2, ndims, 2> flux;

  for (int d = 0; d < ndims; d++) {
    for (int m = 0; m < 2; m++) {
      if (d == 0) {
        flux(d, m) = U(d, k + ks, j + js, i + is + m, n);
      }
      if (d == 1) {
        flux(d, m) = U(d, k + ks, j + js + m, i + is, n);
      }
    }
  }

  // we need to load recon
  if constexpr (std::is_same_v<R, real5d> || std::is_same_v<R, real3d>) {
    SArray<real, 3, ndofs, ndims, 2> recon;
    for (int d = 0; d < ndims; d++) {
      for (int l = 0; l < ndofs; l++) {
        for (int m = 0; m < 2; m++) {
          // full state recon
          if constexpr (std::is_same_v<R, real5d>) {
            if (d == 0) {
              recon(l, d, m) =
                  reconvar(d + l * ndims, k + ks, j + js, i + is + m, n);
            }
            if (d == 1) {
              recon(l, d, m) =
                  reconvar(d + l * ndims, k + ks, j + js + m, i + is, n);
            }
          }
          // reference state
          if constexpr (std::is_same_v<R, real3d>) {
            recon(l, d, m) = reconvar(l, k + ks, n);
          }
        }
      }
    }
    wDnm1bar<ndofs>(tend, recon, flux);
  } else {
    wDnm1bar<ndofs>(tend, reconvar, flux);
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
YAKL_INLINE void compute_wDnm1bar(const real5d &tendvar, const R &reconvar,
                                  const real5d &U, int is, int js, int ks,
                                  int i, int j, int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_wDnm1bar<ndofs>(tend, reconvar, U, is, js, ks, i, j, k, n);

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE Dnm1bar_vert(SArray<real, 1, ndofs> &var,
                              SArray<real, 1, 2> const &flux) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = flux(1) - flux(0);
  }
}

template <uint ndofs>
YAKL_INLINE void compute_Dnm1bar_vert(SArray<real, 1, ndofs> &tend,
                                      const real5d &UW, int is, int js, int ks,
                                      int i, int j, int k, int n) {

  SArray<real, 1, 2> flux;
  for (int m = 0; m < 2; m++) {
    flux(m) = UW(0, k + ks + m, j + js, i + is, n);
  }
  Dnm1bar_vert<ndofs>(tend, flux);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
YAKL_INLINE void compute_Dnm1bar_vert(const real5d &tendvar, const real5d &UW,
                                      int is, int js, int ks, int i, int j,
                                      int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_Dnm1bar_vert<ndofs>(tend, UW, is, js, ks, i, j, k, n);

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE wDnm1bar_vert(SArray<real, 1, ndofs> &var,
                               SArray<real, 2, ndofs, 2> const &recon,
                               SArray<real, 1, 2> const &flux) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = flux(1) * recon(l, 1) - flux(0) * recon(l, 0);
  }
}

template <uint ndofs, class R>
YAKL_INLINE void compute_wDnm1bar_vert(SArray<real, 1, ndofs> &tend,
                                       const R &vertreconvar, const real5d &UW,
                                       int is, int js, int ks, int i, int j,
                                       int k, int n) {

  SArray<real, 1, 2> flux;
  for (int m = 0; m < 2; m++) {
    flux(m) = UW(0, k + ks + m, j + js, i + is, n);
  }

  if constexpr (std::is_same_v<R, real3d> || std::is_same_v<R, real5d>) {
    SArray<real, 2, ndofs, 2> recon;
    for (int m = 0; m < 2; m++) {
      for (int l = 0; l < ndofs; l++) {
        if constexpr (std::is_same_v<R, real5d>) {
          recon(l, m) = vertreconvar(l, k + ks + m, j + js, i + is, n);
        }

        if constexpr (std::is_same_v<R, real3d>) {
          recon(l, m) = vertreconvar(l, k + ks + m, n);
        }
      }
    }
    wDnm1bar_vert<ndofs>(tend, recon, flux);
  } else {
    wDnm1bar_vert<ndofs>(tend, vertreconvar, flux);
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
YAKL_INLINE void compute_wDnm1bar_vert(const real5d &tendvar,
                                       const R &vertreconvar, const real5d &UW,
                                       int is, int js, int ks, int i, int j,
                                       int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_wDnm1bar_vert<ndofs>(tend, vertreconvar, UW, is, js, ks, i, j, k, n);

  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE D1(SArray<real, 1, ndofs> &var,
                    SArray<real, 2, ndofs, 4> const &flux) {

  for (int l = 0; l < ndofs; l++) {

    var(l) = (flux(l, 0) - flux(l, 1) - flux(l, 2) + flux(l, 3));
  }
}

template <uint ndofs>
void YAKL_INLINE compute_D1(SArray<real, 1, ndofs> &tend, const real5d &fluxvar,
                            int is, int js, int ks, int i, int j, int k,
                            int n) {
  SArray<real, 2, ndofs, 4> flux;

  for (int l = 0; l < ndofs; l++) {
    flux(l, 0) = fluxvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    flux(l, 1) = fluxvar(l + 0 * ndofs, k + ks, j + js, i + is, n);
    flux(l, 2) = fluxvar(l + 1 * ndofs, k + ks, j + js, i + is - 1, n);
    flux(l, 3) = fluxvar(l + 0 * ndofs, k + ks, j + js - 1, i + is, n);
  }

  D1<ndofs>(tend, flux);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D1(const real5d &tendvar, const real5d &fluxvar,
                            int is, int js, int ks, int i, int j, int k,
                            int n) {
  SArray<real, 1, ndofs> tend;
  compute_D1<ndofs>(tend, fluxvar, is, js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
  }
  if (addmode == ADD_MODE::ADD) {
    for (int l = 0; l < ndofs; l++) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}

template <uint ndofs>
void YAKL_INLINE compute_D1_ext(SArray<real, 1, ndofs> &tend, const real5d &v,
                                const real5d &w, int is, int js, int ks, int i,
                                int j, int k, int n) {
  SArray<real, 1, 4> flux;
  for (int l = 0; l < ndofs; l++) {
    flux(0) = v(l, k + ks, j + js, i + is, n);     // v1 +
    flux(1) = v(l, k + ks + 1, j + js, i + is, n); // v1 -
    flux(2) = w(l, k + ks, j + js, i + is, n);     // w +
    flux(3) = w(l, k + ks, j + js, i + is - 1, n); // w -
    tend(l) = (flux(0) - flux(1) + flux(2) - flux(3));
  }
}
template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D1_ext(const real5d &tendvar, const real5d &v,
                                const real5d &w, int is, int js, int ks, int i,
                                int j, int k, int n) {
  SArray<real, 1, ndofs> tend;
  compute_D1_ext<ndofs>(tend, v, w, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}

// Fourier transform stuff

void YAKL_INLINE fourier_cwD0(const SArray<complex, 1, ndims> &D0hat,
                              const real c, int i, int j, int k, int nx, int ny,
                              int nz) {
  for (int d = 0; d < ndims; d++) {
    real fac;
    if (d == 0) {
      fac = (2 * pi * i) / nx;
    }
    if (d == 1) {
      fac = (2 * pi * j) / ny;
    }
    // if (d==2) { fac = (2 * pi * k) / nz; }
    complex im(0, 1);
    D0hat(d) = 1._fp - exp(-im * fac);
  }
}

complex YAKL_INLINE fourier_Dnm1bar(const real c, int i, int j, int k, int nx,
                                    int ny, int nz) {
  complex Dnm1barhat;
  for (int d = 0; d < ndims; d++) {
    real fac;
    if (d == 0) {
      fac = (2 * pi * i) / nx;
    }
    if (d == 1) {
      fac = (2 * pi * j) / ny;
    }

    complex im(0, 1);
    Dnm1barhat += exp(im * fac) - 1._fp;
  }
  return Dnm1barhat;
}

void YAKL_INLINE fourier_cwD0Dnm1bar(const SArray<real, 1, ndims> &D0Dnm1barhat,
                                     const real c, int i, int j, int k, int nx,
                                     int ny, int nz) {
  for (int d = 0; d < ndims; d++) {
    real fac;
    if (d == 0) {
      fac = (2 * pi * i) / nx;
    }
    if (d == 1) {
      fac = (2 * pi * j) / ny;
    }
    // if (d==2) { fac = (2 * pi * k) / nz; }
    D0Dnm1barhat(d) = 2 * c * (cos(fac) - 1);
  }
}
