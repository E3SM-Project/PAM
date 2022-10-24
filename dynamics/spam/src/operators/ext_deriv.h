#pragma once

#include "common.h"

complex YAKL_INLINE fourier_Dbar2(const real c, int i, int j, int k, int nx,
                                  int ny, int nz) {
  complex Dbar2hat;
  for (int d = 0; d < ndims; d++) {
    real fac;
    if (d == 0) {
      fac = (2 * pi * i) / nx;
    }
    if (d == 1) {
      fac = (2 * pi * j) / ny;
    }

    complex im(0, 1);
    Dbar2hat += exp(im * fac) - 1._fp;
  }
  return Dbar2hat;
}

void YAKL_INLINE fourier_cwD1Dbar2(const SArray<real, 1, ndims> &D1Dbar2hat,
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
    D1Dbar2hat(d) = 2 * c * (cos(fac) - 1);
  }
}

template <uint ndofs>
void YAKL_INLINE wDbar2(SArray<real, 1, ndofs> &var,
                        SArray<real, 3, ndofs, ndims, 2> const &recon,
                        SArray<real, 2, ndims, 2> const &flux) {

  for (int l = 0; l < ndofs; l++) {
    var(l) = 0.;
    for (int d = 0; d < ndims; d++) {
      var(l) += flux(d, 1) * recon(l, d, 1) - flux(d, 0) * recon(l, d, 0);
    }
  }
}

template <uint ndofs, bool fct, class R>
YAKL_INLINE void compute_wDbar2(SArray<real, 1, ndofs> &tend, const R &reconvar,
                                const real5d &Phivar, const real5d &U, int is,
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
                  reconvar(l + d * ndofs, k + ks, j + js, i + is + m, n);
              if constexpr (fct) {
                recon(l, d, m) *=
                    Phivar(l + d * ndofs, k + ks, j + js, i + is + m, n);
              }
            }
            if (d == 1) {
              recon(l, d, m) =
                  reconvar(l + d * ndofs, k + ks, j + js + m, i + is, n);
              if constexpr (fct) {
                recon(l, d, m) *=
                    Phivar(l + d * ndofs, k + ks, j + js + m, i + is, n);
              }
            }
          }
          // reference state
          if constexpr (std::is_same_v<R, real3d>) {
            recon(l, d, m) = reconvar(l, k, n);
          }
        }
      }
    }
    wDbar2<ndofs>(tend, recon, flux);
  } else {
    wDbar2<ndofs>(tend, reconvar, flux);
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
YAKL_INLINE void compute_wDbar2(const real5d &tendvar, const R &reconvar,
                                const real5d &U, int is, int js, int ks, int i,
                                int j, int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_wDbar2<ndofs, false>(tend, reconvar, U, U, is, js, ks, i, j, k, n);

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

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
YAKL_INLINE void compute_wDbar2_fct(const real5d &tendvar, const R &reconvar,
                                    const real5d &phivar, const real5d &U,
                                    int is, int js, int ks, int i, int j, int k,
                                    int n) {
  SArray<real, 1, ndofs> tend;

  compute_wDbar2<ndofs, true>(tend, reconvar, phivar, U, is, js, ks, i, j, k,
                              n);

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
void YAKL_INLINE wDvbar(SArray<real, 1, ndofs> &var,
                        SArray<real, 2, ndofs, 2> const &recon,
                        SArray<real, 1, 2> const &flux) {
  for (int l = 0; l < ndofs; l++) {
    var(l) = flux(1) * recon(l, 1) - flux(0) * recon(l, 0);
  }
}

template <uint ndofs, bool fct, class R>
YAKL_INLINE void compute_wDvbar(SArray<real, 1, ndofs> &tend,
                                const R &vertreconvar, const real5d &vertphivar,
                                const real5d &UW, int is, int js, int ks, int i,
                                int j, int k, int n) {

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
          recon(l, m) = vertreconvar(l, k + m, n);
        }

        if constexpr (fct) {
          recon(l, m) *= vertphivar(l, k + ks + m, j + js, i + is, n);
        }
      }
    }
    wDvbar<ndofs>(tend, recon, flux);
  } else {
    wDvbar<ndofs>(tend, vertreconvar, flux);
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
YAKL_INLINE void compute_wDvbar(const real5d &tendvar, const R &vertreconvar,
                                const real5d &UW, int is, int js, int ks, int i,
                                int j, int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_wDvbar<ndofs, false>(tend, vertreconvar, UW, UW, is, js, ks, i, j, k,
                               n);

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
template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
YAKL_INLINE void
compute_wDvbar_fct(const real5d &tendvar, const real5d &vertreconvar,
                   const real5d &vertphivar, const real5d &UW, int is, int js,
                   int ks, int i, int j, int k, int n) {
  SArray<real, 1, ndofs> tend;

  compute_wDvbar<ndofs, true>(tend, vertreconvar, vertphivar, UW, is, js, ks, i,
                              j, k, n);

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

void YAKL_INLINE fourier_cwD1(const SArray<complex, 1, ndims> &D1hat,
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
    D1hat(d) = 1._fp - exp(-im * fac);
  }
}

template <uint ndofs>
void YAKL_INLINE wD1(SArray<real, 1, ndims> &var,
                     SArray<real, 2, ndofs, ndims> const &recon,
                     SArray<real, 3, ndofs, ndims, 2> const &dens) {

  for (int d = 0; d < ndims; d++) {
    var(d) = 0.0;
    for (int l = 0; l < ndofs; l++) {
      var(d) += recon(l, d) * (dens(l, d, 1) - dens(l, d, 0));
    }
  }
}

template <uint ndofs, bool fct, class R>
void YAKL_INLINE compute_wD1(SArray<real, 1, ndims> &tend, const R &reconvar,
                             const real5d &Phivar, const real5d &densvar,
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

  // we need to load recon
  if constexpr (std::is_same_v<R, real5d> || std::is_same_v<R, real3d>) {
    SArray<real, 2, ndofs, ndims> recon;
    for (int l = 0; l < ndofs; l++) {
      for (int d = 0; d < ndims; d++) {
        // full state recon
        if constexpr (std::is_same_v<R, real5d>) {
          recon(l, d) = reconvar(l + d * ndofs, k + ks, j + js, i + is, n);
        }
        // reference state
        if constexpr (std::is_same_v<R, real3d>) {
          recon(l, d) = reconvar(l + d * ndofs, k, n);
        }

        // fct
        if constexpr (fct) {
          recon(l, d) *= Phivar(l + d * ndofs, k + ks, j + js, i + is, n);
        }
      }
    }
    wD1<ndofs>(tend, recon, dens);
  } else {
    wD1<ndofs>(tend, reconvar, dens);
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
void YAKL_INLINE compute_wD1(const real5d &tendvar, const R &reconvar,
                             const real5d &densvar, int is, int js, int ks,
                             int i, int j, int k, int n) {
  SArray<real, 1, ndims> tend;
  compute_wD1<ndofs, false>(tend, reconvar, densvar, densvar, is, js, ks, i, j,
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

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_wD1_fct(const real5d &tendvar, const real5d &reconvar,
                                 const real5d &Phivar, const real5d &densvar,
                                 int is, int js, int ks, int i, int j, int k,
                                 int n) {
  SArray<real, 1, ndims> tend;
  compute_wD1<ndofs, true>(tend, Phivar, reconvar, densvar, is, js, ks, i, j, k,
                           n);
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
void YAKL_INLINE wDv(real &var, SArray<real, 1, ndofs> const &recon,
                     SArray<real, 2, ndofs, 2> const &dens) {

  var = 0.0;
  for (int l = 0; l < ndofs; l++) {
    var += recon(l) * (dens(l, 1) - dens(l, 0));
  }
}

template <uint ndofs, bool fct, class R, class D>
real YAKL_INLINE compute_wDv(const R &vertreconvar, const real5d &Phivertvar,
                             const D &densvar, int is, int js, int ks, int i,
                             int j, int k, int n) {
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
      dens(l, 1) = densvar(l, k + 1, n);
      dens(l, 0) = densvar(l, k, n);
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
        recon(l) = vertreconvar(l, k + 1, n);
      }

      if constexpr (fct) {
        recon(l) *= Phivertvar(l, k + ks + 1, j + js, i + is, n);
      }
    }
    wDv<ndofs>(tend, recon, dens);
  } else {
    wDv<ndofs>(tend, vertreconvar, dens);
  }
  return tend;
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R, class D>
void YAKL_INLINE compute_wDv(const real5d &tendvar, const R &vertreconvar,
                             const D &densvar, int is, int js, int ks, int i,
                             int j, int k, int n) {
  real tend = compute_wDv<ndofs, false>(vertreconvar, tendvar, densvar, is, js,
                                        ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    tendvar(0, k + ks, j + js, i + is, n) = tend;
  }
  if (addmode == ADD_MODE::ADD) {
    tendvar(0, k + ks, j + js, i + is, n) += tend;
  }
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE, class R>
void YAKL_INLINE compute_wDv_fct(const real5d &tendvar, const R &vertreconvar,
                                 const real5d &Phivertvar,
                                 const real5d &densvar, int is, int js, int ks,
                                 int i, int j, int k, int n) {
  real tend = compute_wDv<ndofs, true>(vertreconvar, Phivertvar, densvar, is,
                                       js, ks, i, j, k, n);
  if (addmode == ADD_MODE::REPLACE) {
    tendvar(0, k + ks, j + js, i + is, n) = tend;
  }
  if (addmode == ADD_MODE::ADD) {
    tendvar(0, k + ks, j + js, i + is, n) += tend;
  }
}

template <uint ndofs>
void YAKL_INLINE D2(SArray<real, 1, ndofs> &var,
                    SArray<real, 2, ndofs, 4> const &flux) {

  for (int l = 0; l < ndofs; l++) {

    var(l) = (flux(l, 0) - flux(l, 1) - flux(l, 2) + flux(l, 3));
  }
}

template <uint ndofs>
void YAKL_INLINE compute_D2(SArray<real, 1, ndofs> &tend, const real5d &fluxvar,
                            int is, int js, int ks, int i, int j, int k,
                            int n) {
  SArray<real, 2, ndofs, 4> flux;

  for (int l = 0; l < ndofs; l++) {
    flux(l, 0) = fluxvar(l + 1 * ndofs, k + ks, j + js, i + is, n);
    flux(l, 1) = fluxvar(l + 0 * ndofs, k + ks, j + js, i + is, n);
    flux(l, 2) = fluxvar(l + 1 * ndofs, k + ks, j + js, i + is - 1, n);
    flux(l, 3) = fluxvar(l + 0 * ndofs, k + ks, j + js - 1, i + is, n);
  }

  D2<ndofs>(tend, flux);
}

template <uint ndofs, ADD_MODE addmode = ADD_MODE::REPLACE>
void YAKL_INLINE compute_D2(const real5d &tendvar, const real5d &fluxvar,
                            int is, int js, int ks, int i, int j, int k,
                            int n) {
  SArray<real, 1, ndofs> tend;
  compute_D2<ndofs>(tend, fluxvar, is, js, ks, i, j, k, n);
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
void YAKL_INLINE compute_Dxz(SArray<real, 1, ndofs> &tend, const real5d &v,
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
void YAKL_INLINE compute_Dxz(const real5d &tendvar, const real5d &v,
                             const real5d &w, int is, int js, int ks, int i,
                             int j, int k, int n) {
  SArray<real, 1, ndofs> tend;
  compute_Dxz<ndofs>(tend, v, w, is, js, ks, i, j, k, n);
  for (int l = 0; l < ndofs; l++) {
    if (addmode == ADD_MODE::REPLACE) {
      tendvar(l, k + ks, j + js, i + is, n) = tend(l);
    }
    if (addmode == ADD_MODE::ADD) {
      tendvar(l, k + ks, j + js, i + is, n) += tend(l);
    }
  }
}
