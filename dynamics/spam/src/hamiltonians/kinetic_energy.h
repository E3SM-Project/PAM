#pragma once

#include "common.h"
#include "hodge_star.h"
#include "variableset.h"
#include "wedge.h"

namespace pamc {

class Hamiltonian_Hk {

public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  VariableSet varset;
  bool is_initialized;

  Hamiltonian_Hk() { this->is_initialized = false; }

  void initialize(VariableSet &variableset,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->is_initialized = true;
    this->varset = variableset;
  }

  real YAKL_INLINE compute_KE(const real5d &v, const real5d &dens, int is,
                              int js, int ks, int i, int j, int k,
                              int n) const {

    real KE;
    SArray<real, 1, 1> h0, h0im1, h0jm1, h0km1;
    SArray<real, 1, ndims> U, he;
    SArray<real, 2, ndims, 2> h0arr;

    // compute U = H1 v
    compute_H1<1, diff_ord>(U, v, this->primal_geometry, this->dual_geometry,
                            is, js, ks, i, j, k, n);

    // Compute h0 = I h needed for phi calcs
    compute_H2bar<1, diff_ord>(h0, dens, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
    compute_H2bar<1, diff_ord>(h0im1, dens, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i - 1, j, k, n);
    if (ndims >= 2) {
      compute_H2bar<1, diff_ord>(h0jm1, dens, this->primal_geometry,
                                 this->dual_geometry, is, js, ks, i, j - 1, k,
                                 n);
    }
    if (ndims >= 3) {
      compute_H2bar<1, diff_ord>(h0km1, dens, this->primal_geometry,
                                 this->dual_geometry, is, js, ks, i, j, k - 1,
                                 n);
    }

    // compute he = phi h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        h0arr(d, 0) = h0(0);
        h0arr(d, 1) = h0im1(0);
      }
      if (d == 1) {
        h0arr(d, 0) = h0(0);
        h0arr(d, 1) = h0jm1(0);
      }
      if (d == 2) {
        h0arr(d, 0) = h0(0);
        h0arr(d, 1) = h0km1(0);
      }
    }
    phi(he, h0arr);

    KE = 0.;
    for (int d = 0; d < ndims; d++) {
      KE = KE + he(d) * (U(d) * v(d, k + ks, j + js, i + is, n));
    }
    return 0.5_fp * KE;
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_F_and_K(const real5d &F, const real5d &K,
                                   const real5d &v, const real5d &U,
                                   const real5d &dens0, int is, int js, int ks,
                                   int i, int j, int k, int n,
                                   real fac = 1._fp) const {
    SArray<real, 2, ndims, 2> D0;
    SArray<real, 1, ndims> he;

    // compute he = phi * h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js - 1, i + is, n);
      }
      if (d == 2) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks - 1, j + js, i + is, n);
      }
    }
    phi(he, D0);

    // compute F = he * U
    for (int d = 0; d < ndims; d++) {
      if (addmode == ADD_MODE::REPLACE) {
        F(d, k + ks, j + js, i + is, n) =
            fac * U(d, k + ks, j + js, i + is, n) * he(d);
      } else if (addmode == ADD_MODE::ADD) {
        F(d, k + ks, j + js, i + is, n) +=
            fac * U(d, k + ks, j + js, i + is, n) * he(d);
      }
    }

    // compute K = 1/2 * PhiT(U,V)
    compute_phiT(K, U, v, is, js, ks, i, j, k, n);
    K(0, k + ks, j + js, i + is, n) *= 0.5_fp;
  }

  void YAKL_INLINE compute_F_and_he(const real5d &F, const real5d &HE,
                                    const real5d &U, const real5d &dens0,
                                    int is, int js, int ks, int i, int j, int k,
                                    int n) const {
    SArray<real, 2, ndims, 2> D0;
    SArray<real, 1, ndims> he;

    // compute he = phi * h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js - 1, i + is, n);
      }
      if (d == 2) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks - 1, j + js, i + is, n);
      }
    }
    phi(he, D0);

    // compute F = he * U, set HE
    for (int d = 0; d < ndims; d++) {
      F(d, k + ks, j + js, i + is, n) = U(d, k + ks, j + js, i + is, n) * he(d);
      HE(d, k + ks, j + js, i + is, n) = he(d);
    }
  }

  void YAKL_INLINE compute_he(const real5d &HE, const real5d &dens0, int is,
                              int js, int ks, int i, int j, int k,
                              int n) const {
    SArray<real, 2, ndims, 2> D0;
    SArray<real, 1, ndims> he;

    // compute he = phi * h0
    for (int d = 0; d < ndims; d++) {
      if (d == 0) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js, i + is - 1, n);
      }
      if (d == 1) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks, j + js - 1, i + is, n);
      }
      if (d == 2) {
        D0(d, 0) = dens0(0, k + ks, j + js, i + is, n);
        D0(d, 1) = dens0(0, k + ks - 1, j + js, i + is, n);
      }
    }
    phi(he, D0);

    // set HE
    for (int d = 0; d < ndims; d++) {
      HE(d, k + ks, j + js, i + is, n) = he(d);
    }
  }

  // FIX THIS TO GET TOTAL DENSITY FROM VARSET!
  //  Note that this ADDS to Bvar...
  void YAKL_INLINE compute_dKddens(const real5d &B, const real5d &K, int is,
                                   int js, int ks, int i, int j, int k, int n,
                                   real fac = 1._fp) const {
    SArray<real, 1, 1> K0;
    compute_H2bar<1, diff_ord>(K0, K, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
    B(varset.active_id_mass, k + ks, j + js, i + is, n) += fac * K0(0);
  }
};

class Hamiltonian_Hk_extruded {

public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  VariableSet varset;
  bool is_initialized;

  Hamiltonian_Hk_extruded() { this->is_initialized = false; }

  void initialize(VariableSet &variableset,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->is_initialized = true;
    this->varset = variableset;
  }

  real YAKL_INLINE compute_KE_top(const real5d &v, const real5d &w,
                                  const real5d &dens, int is, int js, int ks,
                                  int i, int j, int k, int n) const {
    real K2 = 0.;
    SArray<real, 1, 1> UW0;
    compute_H01<1, vert_diff_ord>(UW0, w, this->primal_geometry,
                                  this->dual_geometry, is, js, ks, i, j, k, n);
    real w0;
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w0 = w(0, k + ks - 1, j + js, i + is, n);
    K2 += 0.5 * w0 * UW0(0);
    return _compute_KE(K2, v, w, dens, is, js, ks, i, j, k, n);
  }

  real YAKL_INLINE compute_KE_bottom(const real5d &v, const real5d &w,
                                     const real5d &dens, int is, int js, int ks,
                                     int i, int j, int k, int n) const {
    real K2 = 0.;
    SArray<real, 1, 1> UW1;
    compute_H01<1, vert_diff_ord>(UW1, w, this->primal_geometry,
                                  this->dual_geometry, is, js, ks, i, j, k + 1,
                                  n);
    real w1;
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w1 = w(0, k + ks, j + js, i + is, n);
    K2 += 0.5 * w1 * UW1(0);
    return _compute_KE(K2, v, w, dens, is, js, ks, i, j, k, n);
  }

  real YAKL_INLINE compute_KE(const real5d &v, const real5d &w,
                              const real5d &dens, int is, int js, int ks, int i,
                              int j, int k, int n) const {
    real K2 = 0.;
    SArray<real, 1, 1> UW0, UW1;
    compute_H01<1, vert_diff_ord>(UW0, w, this->primal_geometry,
                                  this->dual_geometry, is, js, ks, i, j, k, n);
    compute_H01<1, vert_diff_ord>(UW1, w, this->primal_geometry,
                                  this->dual_geometry, is, js, ks, i, j, k + 1,
                                  n);
    real w0, w1;
    // Have to subtract 1 from k here since UW has an extra dof compared to w
    w0 = w(0, k + ks - 1, j + js, i + is, n);
    w1 = w(0, k + ks, j + js, i + is, n);
    K2 += 0.5 * (w0 * UW0(0) + w1 * UW1(0));
    return _compute_KE(K2, v, w, dens, is, js, ks, i, j, k, n);
  }

  real YAKL_INLINE _compute_KE(real K2, const real5d &v, const real5d &w,
                               const real5d &dens, int is, int js, int ks,
                               int i, int j, int k, int n) const {
    real v0, v1;
    SArray<real, 1, ndims> U0, U1;
    compute_H10<1, diff_ord>(U0, v, this->primal_geometry, this->dual_geometry,
                             is, js, ks, i, j, k, n);
    compute_H10<1, diff_ord>(U1, v, this->primal_geometry, this->dual_geometry,
                             is, js, ks, i + 1, j, k, n);
    v0 = v(0, k + ks, j + js, i + is, n);
    v1 = v(0, k + ks, j + js, i + is + 1, n);
    K2 += 0.5 * (v0 * U0(0) + v1 * U1(0));

    if (ndims == 2) {
      real v0, v1;
      SArray<real, 1, ndims> U0, U1;
      compute_H10<1, diff_ord>(U0, v, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
      compute_H10<1, diff_ord>(U1, v, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j + 1, k, n);
      v0 = v(1, k + ks, j + js, i + is, n);
      v1 = v(1, k + ks, j + js + 1, i + is, n);
      K2 += 0.5 * (v0 * U0(1) + v1 * U1(1));
    }

    K2 *= 0.5;

    // Compute h0 = I h needed for phi calcs
    SArray<real, 1, 1> h0;

    const auto total_density_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      return varset.get_total_density(densvar, k, j, i, 0, 0, 0, n);
    };
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        total_density_f, h0, dens, this->primal_geometry, this->dual_geometry,
        is, js, ks, i, j, k, n);

    return h0(0) * K2;
  }

  YAKL_INLINE void compute_he_U_and_K(SArray<real, 1, ndims> &he, real &hew,
                                      SArray<real, 1, ndims> &U, real &UW,
                                      real &K2, real5d densvar, real5d Vvar,
                                      real5d Wvar, int is, int js, int ks,
                                      int i, int j, int k, int n) const {

#if defined PAMC_AN || defined PAMC_MAN
    auto &rho_pi = varset.reference_state.rho_pi.data;
    auto &rho_di = varset.reference_state.rho_di.data;
    for (int d = 0; d < ndims; ++d) {
      he(d) = rho_pi(0, k + ks, n);
    }
    hew = rho_di(0, k + ks, n);
#else
    SArray<real, 1, 1> dens0_ijk, dens0_im1, dens0_jm1, dens0_km1;

    const auto total_density_f =
        YAKL_LAMBDA(const real5d &densvar, int d, int k, int j, int i, int n) {
      return varset.get_total_density(densvar, k, j, i, 0, 0, 0, n);
    };

    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        total_density_f, dens0_ijk, densvar, primal_geometry, dual_geometry, is,
        js, ks, i, j, k, n);
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        total_density_f, dens0_im1, densvar, primal_geometry, dual_geometry, is,
        js, ks, i - 1, j, k, n);
    if (ndims > 1) {
      compute_Hn1bar<1, diff_ord, vert_diff_ord>(
          total_density_f, dens0_jm1, densvar, primal_geometry, dual_geometry,
          is, js, ks, i, j - 1, k, n);
    }
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        total_density_f, dens0_km1, densvar, primal_geometry, dual_geometry, is,
        js, ks, i, j, k - 1, n);

    he(0) = 0.5_fp * (dens0_ijk(0) + dens0_im1(0));
    if (ndims > 1) {
      he(1) = 0.5_fp * (dens0_ijk(0) + dens0_jm1(0));
    }
    hew = 0.5_fp * (dens0_ijk(0) + dens0_km1(0));
#endif

    SArray<real, 1, ndims> u_ijk, u_ip1, u_jp1;
    SArray<real, 1, 1> uw_ijk, uw_kp1;

    compute_H10<1, diff_ord>(u_ijk, Vvar, primal_geometry, dual_geometry, is,
                             js, ks, i, j, k, n);
    compute_H10<1, diff_ord>(u_ip1, Vvar, primal_geometry, dual_geometry, is,
                             js, ks, i + 1, j, k, n);
    if (ndims > 1) {
      compute_H10<1, diff_ord>(u_jp1, Vvar, primal_geometry, dual_geometry, is,
                               js, ks, i, j + 1, k, n);
    }

    const int dni = dual_geometry.topology.ni;
    const int dnl = dual_geometry.topology.nl;

    if (k == 0 || k == (dni - 1)) {
      uw_ijk(0) = 0;
    } else {
      compute_H01<1, vert_diff_ord>(uw_ijk, Wvar, primal_geometry,
                                    dual_geometry, is, js, ks, i, j, k, n);
    }

    if (k >= (dni - 2)) {
      uw_kp1(0) = 0;
    } else {
      compute_H01<1, vert_diff_ord>(uw_kp1, Wvar, primal_geometry,
                                    dual_geometry, is, js, ks, i, j, k + 1, n);
    }

    for (int d = 0; d < ndims; ++d) {
      U(d) = u_ijk(d);
    }
    UW = uw_ijk(0);

    K2 = 0.5_fp * (Vvar(0, k + ks, j + js, i + is, n) * u_ijk(0) +
                   Vvar(0, k + ks, j + js, i + 1 + is, n) * u_ip1(0));
    if (ndims > 1) {
      K2 += 0.5_fp * (Vvar(1, k + ks, j + js, i + is, n) * u_ijk(1) +
                      Vvar(1, k + ks, j + 1 + js, i + is, n) * u_jp1(1));
    }
    if (k < dnl) {
      K2 += 0.5_fp * (Wvar(0, k - 1 + ks, j + js, i + is, n) * uw_ijk(0) +
                      Wvar(0, k + ks, j + js, i + is, n) * uw_kp1(0));
    }

    K2 *= 0.5;
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dKddens(const real5d &B, const real5d &K, int is,
                                   int js, int ks, int i, int j, int k, int n,
                                   real fac = 1._fp) const {
    SArray<real, 1, 1> K0;
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(K0, K, this->primal_geometry,
                                               this->dual_geometry, is, js, ks,
                                               i, j, k, n);
    if (addmode == ADD_MODE::REPLACE) {
      B(varset.active_id_mass, k + ks, j + js, i + is, n) = fac * K0(0);
    }
    if (addmode == ADD_MODE::ADD) {
      B(varset.active_id_mass, k + ks, j + js, i + is, n) += fac * K0(0);
    }
  }
};
} // namespace pamc
