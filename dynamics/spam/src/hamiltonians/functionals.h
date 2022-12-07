
#pragma once

#include "common.h"
#include "ext_deriv.h"
#include "variableset.h"
#include "wedge.h"

struct pvpe {
  real pv = 0., pe = 0.;
};

// ADD SUPPORT FOR GETTING TOTAL DENSITY FROM VARSET...
class Functional_PVPE {
public:
  bool is_initialized;
  VariableSet varset;

  Functional_PVPE() { this->is_initialized = false; }

  void initialize(VariableSet &variableset) {
    this->is_initialized = true;
    this->varset = variableset;
  }

  real YAKL_INLINE compute_hv(const real5d &dens, int is, int js, int ks, int i,
                              int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    Dv(0) = dens(0, k + ks, j + js, i + is, n);
    Dv(1) = dens(0, k + ks, j + js, i + is - 1, n);
    Dv(2) = dens(0, k + ks, j + js - 1, i + is, n);
    Dv(3) = dens(0, k + ks, j + js - 1, i + is - 1, n);
    R(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_zeta(const real5d &v, int is, int js, int ks, int i,
                                int j, int k, int n) const {
    SArray<real, 1, 1> zeta;
    // compute zeta = D1 v
    compute_D1<1>(zeta, v, is, js, ks, i, j, k, n);
    return zeta(0);
  }

  real YAKL_INLINE compute_eta(const real5d &v, const real5d &coriolis, int is,
                               int js, int ks, int i, int j, int k,
                               int n) const {
    real zeta = compute_zeta(v, is, js, ks, i, j, k, n);
    return zeta + coriolis(0, k + ks, j + js, i + is, n);
  }

  // This computes relative q0
  void YAKL_INLINE compute_q0f0(const real5d &q0, const real5d &f0,
                                const real5d &v, const real5d &dens,
                                const real5d &coriolis, int is, int js, int ks,
                                int i, int j, int k, int n) const {

    real hv = compute_hv(dens, is, js, ks, i, j, k, n);
    real zeta = compute_zeta(v, is, js, ks, i, j, k, n);

    // compute q0 = zeta / hv and f0 = f / hv
    q0(0, k + ks, j + js, i + is, n) = zeta / hv;
    f0(0, k + ks, j + js, i + is, n) =
        coriolis(0, k + ks, j + js, i + is, n) / hv;
  }

  // This computes TRUE q0
  void YAKL_INLINE compute_q0(const real5d &q0, const real5d &v,
                              const real5d &dens, const real5d &coriolis,
                              int is, int js, int ks, int i, int j, int k,
                              int n) const {
    real hv = compute_hv(dens, is, js, ks, i, j, k, n);
    real eta = compute_eta(v, coriolis, is, js, ks, i, j, k, n);
    // compute q0 = zeta / hv and f0 = f / hv
    q0(0, k + ks, j + js, i + is, n) = eta / hv;
  }

  pvpe YAKL_INLINE compute_PVPE(const real5d &v, const real5d &dens,
                                const real5d &coriolis, int is, int js, int ks,
                                int i, int j, int k, int n) const {
    pvpe vals;
    real eta = compute_eta(v, coriolis, is, js, ks, i, j, k, n);
    real hv = compute_hv(dens, is, js, ks, i, j, k, n);
    real q0 = eta / hv;

    vals.pv = eta;
    vals.pe = 0.5_fp * eta * q0;

    return vals;
  }
};

//
// class Functional_PVPE_rhod {
// public:
//   bool is_initialized;
//
//    Functional_PVPE_rhod() {
//      this->is_initialized = false;
// }
//
// void initialize(Parameters &params)
// {
//   this->is_initialized = true;
// }
//
// real YAKL_INLINE compute_hv(const real5d& dens, const real5d& densfct, int
// is, int js, int ks, int i, int j, int k)
// {
//   SArray<real,1> hv;
//   SArray<real,4> Dv;
//
//   // compute hv = R h
//   // Uses linearity of R
//   Dv(0) = dens(0, k+ks, j+js, i+is) + densfct(0, k+ks, j+js, i+is) +
//   densfct(1, k+ks, j+js, i+is) + densfct(2, k+ks, j+js, i+is); Dv(1) =
//   dens(0, k+ks, j+js, i+is-1) + densfct(0, k+ks, j+js, i+is-1) + densfct(1,
//   k+ks, j+js, i+is-1) + densfct(2, k+ks, j+js, i+is-1); Dv(2) = dens(0, k+ks,
//   j+js-1, i+is) + densfct(0, k+ks, j+js-1, i+is) + densfct(1, k+ks, j+js-1,
//   i+is) + densfct(2, k+ks, j+js-1, i+is); Dv(3) = dens(0, k+ks, j+js-1,
//   i+is-1) + densfct(0, k+ks, j+js-1, i+is-1) + densfct(1, k+ks, j+js-1,
//   i+is-1) + densfct(2, k+ks, j+js-1, i+is-1); R(hv, Dv);
//
//   return hv(0);
// }
//
// real YAKL_INLINE compute_zeta(const real5d& v, int is, int js, int ks, int i,
// int j, int k)
// {
//   SArray<real,1> zeta;
//   // compute zeta = D2 v
//   compute_D2<1>(zeta, v, is, js, ks, i, j, k);
//   return zeta(0);
// }
//
// real YAKL_INLINE compute_eta(const real5d& v, const real5d& coriolis, int is,
// int js, int ks, int i, int j, int k)
// {
//   real zeta = compute_zeta(v, is, js, ks, i, j, k);
//   return zeta + coriolis(0, k+ks, j+js, i+is);
// }
//
//
// void YAKL_INLINE compute_q0f0(const real5d& q0, const real5d& f0, const
// real5d& v, const const real5d& dens, const real5d& densfct, const real5d&
// coriolis, int is, int js, int ks, int i, int j, int k)
// {
//
//   real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
//   real zeta = compute_zeta(v, is, js, ks, i, j, k);
//
//   // compute q0 = zeta / hv and f0 = f / hv
//     q0(0, k+ks, j+js, i+is) = zeta / hv;
//     f0(0, k+ks, j+js, i+is) = coriolis(0, k+ks, j+js, i+is) / hv;
//
// }
//
// void YAKL_INLINE compute_q0(const real5d& q0, const real5d& v, const real5d&
// dens, const real5d& densfct, int is, int js, int ks, int i, int j, int k)
// {
// real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
// real zeta = compute_zeta(v, is, js, ks, i, j, k);
// // compute q0 = zeta / hv and f0 = f / hv
//   q0(0, k+ks, j+js, i+is) = zeta / hv;
// }
//
// pvpe YAKL_INLINE compute_PVPE(const real5d& v, const real5d& dens, const
// real5d& densfct, const real5d& coriolis, int is, int js, int ks, int i, int
// j, int k)
// {
//   pvpe vals;
//   real eta = compute_eta(v, coriolis, is, js, ks, i, j, k);
//   real hv = compute_hv(dens, densfct, is, js, ks, i, j, k);
//   real q0 = eta / hv;
//
//   vals.pv = eta;
//   vals.pe = 0.5 * eta * q0;
//
//   return vals;
// }
//
// };
//
//

class Functional_PVPE_extruded {

public:
  bool is_initialized;
  VariableSet varset;

  Functional_PVPE_extruded() { this->is_initialized = false; }

  void initialize(VariableSet &variableset) {
    this->is_initialized = true;
    this->varset = variableset;
  }

  real YAKL_INLINE compute_hvxz(const real5d &dens, int is, int js, int ks,
                                int i, int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = dens(0, k + ks, j + js, i + is, n);
    Dv(1) = dens(0, k + ks, j + js, i + is - 1, n);
    Dv(2) = dens(0, k + ks + 1, j + js, i + is, n);
    Dv(3) = dens(0, k + ks + 1, j + js, i + is - 1, n);
    R(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvxz_top(const real5d &dens, int is, int js, int ks,
                                    int i, int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = dens(0, k + ks, j + js, i + is, n);
    Dv(1) = dens(0, k + ks, j + js, i + is - 1, n);
    Dv(2) = dens(0, k + ks + 1, j + js, i + is, n);     // gets 1/2
    Dv(3) = dens(0, k + ks + 1, j + js, i + is - 1, n); // gets 1/2
    Rbnd(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvxz_bottom(const real5d &dens, int is, int js,
                                       int ks, int i, int j, int k,
                                       int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = dens(0, k + ks + 1, j + js, i + is, n);
    Dv(1) = dens(0, k + ks + 1, j + js, i + is - 1, n);
    Dv(2) = dens(0, k + ks, j + js, i + is, n);     // gets 1/2
    Dv(3) = dens(0, k + ks, j + js, i + is - 1, n); // gets 1/2
    Rbnd(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_zetaxz(const real5d &v, const real5d &w, int is,
                                  int js, int ks, int i, int j, int k,
                                  int n) const {
    SArray<real, 1, 1> zeta;
    // compute zeta = D1_ext "v"
    compute_D1_ext<1>(zeta, v, w, is, js, ks, i, j, k, n);
    return zeta(0);
  }

  real YAKL_INLINE compute_etaxz(const real5d &v, const real5d &w,
                                 const real5d &coriolisxz, int is, int js,
                                 int ks, int i, int j, int k, int n) const {
    real zeta = compute_zetaxz(v, w, is, js, ks, i, j, k, n);
    return zeta + coriolisxz(0, k + ks, j + js, i + is, n);
  }

  // This computes true qxz
  void YAKL_INLINE compute_qxz0(const real5d &qxz0, const real5d &v,
                                const real5d &w, const real5d &dens,
                                const real5d &coriolisxz, int is, int js,
                                int ks, int i, int j, int k, int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    real hv = compute_hvxz(dens, is, js, ks, i, j, k - 1, n);
    real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    qxz0(0, k + ks, j + js, i + is, n) = eta / hv;
  }

  // This computes true qxz
  void YAKL_INLINE compute_qxz0_top(const real5d &qxz0, const real5d &v,
                                    const real5d &w, const real5d &dens,
                                    const real5d &coriolisxz, int is, int js,
                                    int ks, int i, int j, int k, int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    real hv = compute_hvxz_top(dens, is, js, ks, i, j, k - 1, n);
    real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    qxz0(0, k + ks, j + js, i + is, n) = eta / hv;
  }

  // This computes true qxz
  void YAKL_INLINE compute_qxz0_bottom(const real5d &qxz0, const real5d &v,
                                       const real5d &w, const real5d &dens,
                                       const real5d &coriolisxz, int is, int js,
                                       int ks, int i, int j, int k,
                                       int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    real hv = compute_hvxz_bottom(dens, is, js, ks, i, j, k - 1, n);
    real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    qxz0(0, k + ks, j + js, i + is, n) = eta / hv;
  }

  // This computes relative qxz
  void YAKL_INLINE compute_qxz0fxz0(const real5d &qxz0, const real5d &fxz0,
                                    const real5d &v, const real5d &w,
                                    const real5d &dens,
                                    const real5d &coriolisxz, int is, int js,
                                    int ks, int i, int j, int k, int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    real hv = compute_hvxz(dens, is, js, ks, i, j, k - 1, n);
    real zeta = compute_zetaxz(v, w, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    qxz0(0, k + ks, j + js, i + is, n) = zeta / hv;
    fxz0(0, k + ks, j + js, i + is, n) =
        coriolisxz(0, k + ks, j + js, i + is, n) / hv;
  }

  // This computes relative qxz
  void YAKL_INLINE compute_qxz0fxz0_top(const real5d &qxz0, const real5d &fxz0,
                                        const real5d &v, const real5d &w,
                                        const real5d &dens,
                                        const real5d &coriolisxz, int is,
                                        int js, int ks, int i, int j, int k,
                                        int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    real hv = compute_hvxz_top(dens, is, js, ks, i, j, k - 1, n);
    real zeta = compute_zetaxz(v, w, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    qxz0(0, k + ks, j + js, i + is, n) = zeta / hv;
    fxz0(0, k + ks, j + js, i + is, n) =
        coriolisxz(0, k + ks, j + js, i + is, n) / hv;
  }

  // This computes relative qxz
  void YAKL_INLINE compute_qxz0fxz0_bottom(const real5d &qxz0,
                                           const real5d &fxz0, const real5d &v,
                                           const real5d &w, const real5d &dens,
                                           const real5d &coriolisxz, int is,
                                           int js, int ks, int i, int j, int k,
                                           int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    real hv = compute_hvxz_bottom(dens, is, js, ks, i, j, k - 1, n);
    real zeta = compute_zetaxz(v, w, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    qxz0(0, k + ks, j + js, i + is, n) = zeta / hv;
    fxz0(0, k + ks, j + js, i + is, n) =
        coriolisxz(0, k + ks, j + js, i + is, n) / hv;
  }

  pvpe YAKL_INLINE compute_PVPE(const real5d &v, const real5d &w,
                                const real5d &dens, const real5d &coriolisxz,
                                int is, int js, int ks, int i, int j, int k,
                                int n) const {
    pvpe vals;
    // No subtraction here since this is called on primal cells p11
    real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k, n);
    real hv = compute_hvxz(dens, is, js, ks, i, j, k, n);
    real q0 = eta / hv;
    vals.pv = eta;
    vals.pe = 0.5_fp * eta * q0;
    return vals;
  }

  pvpe YAKL_INLINE compute_PVPE_top(const real5d &v, const real5d &w,
                                    const real5d &dens,
                                    const real5d &coriolisxz, int is, int js,
                                    int ks, int i, int j, int k, int n) const {
    pvpe vals;
    // No subtraction here since this is called on primal cells p11
    real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k, n);
    real hv = compute_hvxz_top(dens, is, js, ks, i, j, k, n);
    real q0 = eta / hv;
    vals.pv = eta;
    vals.pe = 0.5_fp * eta * q0;
    return vals;
  }

  pvpe YAKL_INLINE compute_PVPE_bottom(const real5d &v, const real5d &w,
                                       const real5d &dens,
                                       const real5d &coriolisxz, int is, int js,
                                       int ks, int i, int j, int k,
                                       int n) const {
    pvpe vals;
    // No subtraction here since this is called on primal cells p11
    real eta = compute_etaxz(v, w, coriolisxz, is, js, ks, i, j, k, n);
    real hv = compute_hvxz_bottom(dens, is, js, ks, i, j, k, n);
    real q0 = eta / hv;
    vals.pv = eta;
    vals.pe = 0.5_fp * eta * q0;
    return vals;
  }
};
