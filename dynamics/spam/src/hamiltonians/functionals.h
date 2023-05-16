
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
    Dv(0) = varset.get_total_density(dens, k, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k, j, i - 1, ks, js, is, n);
    Dv(2) = varset.get_total_density(dens, k + 1, j, i, ks, js, is, n);
    Dv(3) = varset.get_total_density(dens, k + 1, j, i - 1, ks, js, is, n);
    R(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvxz_top(const real5d &dens, int is, int js, int ks,
                                    int i, int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = varset.get_total_density(dens, k, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k, j, i - 1, ks, js, is, n);
    Dv(2) =
        varset.get_total_density(dens, k + 1, j, i, ks, js, is, n); // gets 1/2
    Dv(3) = varset.get_total_density(dens, k + 1, j, i - 1, ks, js, is,
                                     n); // gets 1/2
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
    Dv(0) = varset.get_total_density(dens, k + 1, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k + 1, j, i - 1, ks, js, is, n);
    Dv(2) = varset.get_total_density(dens, k, j, i, ks, js, is, n); // gets 1/2
    Dv(3) =
        varset.get_total_density(dens, k, j, i - 1, ks, js, is, n); // gets 1/2
    Rbnd(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvyz(const real5d &dens, int is, int js, int ks,
                                int i, int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = varset.get_total_density(dens, k, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k, j - 1, i, ks, js, is, n);
    Dv(2) = varset.get_total_density(dens, k + 1, j, i, ks, js, is, n);
    Dv(3) = varset.get_total_density(dens, k + 1, j - 1, i, ks, js, is, n);
    R(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvyz_top(const real5d &dens, int is, int js, int ks,
                                    int i, int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = varset.get_total_density(dens, k, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k, j - 1, i, ks, js, is, n);
    Dv(2) =
        varset.get_total_density(dens, k + 1, j, i, ks, js, is, n); // gets 1/2
    Dv(3) = varset.get_total_density(dens, k + 1, j - 1, i, ks, js, is,
                                     n); // gets 1/2
    Rbnd(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvyz_bottom(const real5d &dens, int is, int js,
                                       int ks, int i, int j, int k,
                                       int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = varset.get_total_density(dens, k + 1, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k + 1, j - 1, i, ks, js, is, n);
    Dv(2) = varset.get_total_density(dens, k, j, i, ks, js, is, n); // gets 1/2
    Dv(3) =
        varset.get_total_density(dens, k, j - 1, i, ks, js, is, n); // gets 1/2
    Rbnd(hv, Dv);

    return hv(0);
  }

  real YAKL_INLINE compute_hvxy(const real5d &dens, int is, int js, int ks,
                                int i, int j, int k, int n) const {
    SArray<real, 1, 1> hv;
    SArray<real, 1, 4> Dv;

    // compute hv = R h
    // Uses linearity of R
    Dv(0) = varset.get_total_density(dens, k, j, i, ks, js, is, n);
    Dv(1) = varset.get_total_density(dens, k, j, i - 1, ks, js, is, n);
    Dv(2) = varset.get_total_density(dens, k, j - 1, i, ks, js, is, n);
    Dv(3) = varset.get_total_density(dens, k, j - 1, i - 1, ks, js, is, n);
    R(hv, Dv);

    return hv(0);
  }

  void YAKL_INLINE compute_zetahz(SArray<real, 1, ndims> &zetahz,
                                  const real5d &v, const real5d &w, int is,
                                  int js, int ks, int i, int j, int k,
                                  int n) const {
    SArray<real, 2, 1, ndims> zeta;
    // compute zeta = D1_ext "v"
    compute_D1_ext<1>(zeta, v, w, is, js, ks, i, j, k, n);
    for (int d = 0; d < ndims; ++d) {
      zetahz(d) = zeta(0, d);
    }
  }

  void YAKL_INLINE compute_etahz(SArray<real, 1, ndims> &etahz, const real5d &v,
                                 const real5d &w, const real5d &coriolishz,
                                 int is, int js, int ks, int i, int j, int k,
                                 int n) const {
    SArray<real, 1, ndims> zetahz;
    compute_zetahz(zetahz, v, w, is, js, ks, i, j, k, n);
    for (int d = 0; d < ndims; ++d) {
      etahz(d) = zetahz(d) + coriolishz(d, k + ks, j + js, i + is, n);
    }
  }

  // This computes true qhz
  void YAKL_INLINE compute_qhz(const real5d &qhz, const real5d &v,
                               const real5d &w, const real5d &dens,
                               const real5d &coriolishz, int is, int js, int ks,
                               int i, int j, int k, int n) const {
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    SArray<real, 1, ndims> hv;
    hv(0) = compute_hvxz(dens, is, js, ks, i, j, k - 1, n);
    if (ndims > 1) {
      hv(1) = compute_hvyz(dens, is, js, ks, i, j, k - 1, n);
    }

    SArray<real, 1, ndims> etahz;
    compute_etahz(etahz, v, w, coriolishz, is, js, ks, i, j, k - 1, n);

    // compute q0 = zeta / hv and f0 = f / hv
    for (int d = 0; d < ndims; ++d) {
      qhz(d, k + ks, j + js, i + is, n) = etahz(d) / hv(d);
    }
  }

  // This computes true qhz
  void YAKL_INLINE compute_qhz_top(const real5d &qhz, const real5d &v,
                                   const real5d &w, const real5d &dens,
                                   const real5d &coriolishz, int is, int js,
                                   int ks, int i, int j, int k, int n) const {
    SArray<real, 1, ndims> hv;
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    hv(0) = compute_hvxz_top(dens, is, js, ks, i, j, k - 1, n);
    if (ndims > 1) {
      hv(1) = compute_hvyz_top(dens, is, js, ks, i, j, k - 1, n);
    }
    SArray<real, 1, ndims> etahz;
    compute_etahz(etahz, v, w, coriolishz, is, js, ks, i, j, k - 1, n);

    // compute q0 = zeta / hv and f0 = f / hv
    for (int d = 0; d < ndims; ++d) {
      qhz(d, k + ks, j + js, i + is, n) = etahz(d) / hv(d);
    }
  }

  // This computes true qhz
  void YAKL_INLINE compute_qhz_bottom(const real5d &qhz, const real5d &v,
                                      const real5d &w, const real5d &dens,
                                      const real5d &coriolishz, int is, int js,
                                      int ks, int i, int j, int k,
                                      int n) const {
    SArray<real, 1, ndims> hv;
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    hv(0) = compute_hvxz_bottom(dens, is, js, ks, i, j, k - 1, n);
    if (ndims > 1) {
      hv(1) = compute_hvyz_bottom(dens, is, js, ks, i, j, k - 1, n);
    }
    SArray<real, 1, ndims> etahz;
    compute_etahz(etahz, v, w, coriolishz, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    for (int d = 0; d < ndims; ++d) {
      qhz(d, k + ks, j + js, i + is, n) = etahz(d) / hv(d);
    }
  }

  // This computes relative qhz
  void YAKL_INLINE compute_qhzfhz(const real5d &qhz, const real5d &fhz,
                                  const real5d &v, const real5d &w,
                                  const real5d &dens, const real5d &coriolishz,
                                  int is, int js, int ks, int i, int j, int k,
                                  int n) const {
    SArray<real, 1, ndims> hv;
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    hv(0) = compute_hvxz(dens, is, js, ks, i, j, k - 1, n);
    if (ndims > 1) {
      hv(1) = compute_hvyz(dens, is, js, ks, i, j, k - 1, n);
    }
    SArray<real, 1, ndims> zetahz;
    compute_zetahz(zetahz, v, w, is, js, ks, i, j, k - 1, n);
    // compute q0 = zeta / hv and f0 = f / hv
    for (int d = 0; d < ndims; ++d) {
      qhz(d, k + ks, j + js, i + is, n) = zetahz(d) / hv(d);
      fhz(d, k + ks, j + js, i + is, n) =
          coriolishz(d, k + ks, j + js, i + is, n) / hv(d);
    }
  }

  // This computes relative qhz
  void YAKL_INLINE compute_qhzfhz_top(const real5d &qhz, const real5d &fhz,
                                      const real5d &v, const real5d &w,
                                      const real5d &dens,
                                      const real5d &coriolishz, int is, int js,
                                      int ks, int i, int j, int k,
                                      int n) const {
    SArray<real, 1, ndims> hv;
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    hv(0) = compute_hvxz_top(dens, is, js, ks, i, j, k - 1, n);
    if (ndims > 1) {
      hv(1) = compute_hvyz_top(dens, is, js, ks, i, j, k - 1, n);
    }

    SArray<real, 1, ndims> zetahz;
    compute_zetahz(zetahz, v, w, is, js, ks, i, j, k - 1, n);

    // compute q0 = zeta / hv and f0 = f / hv
    for (int d = 0; d < ndims; ++d) {
      qhz(d, k + ks, j + js, i + is, n) = zetahz(d) / hv(d);
      fhz(d, k + ks, j + js, i + is, n) =
          coriolishz(d, k + ks, j + js, i + is, n) / hv(d);
    }
  }

  // This computes relative qhz
  void YAKL_INLINE compute_qhzfhz_bottom(const real5d &qhz, const real5d &fhz,
                                         const real5d &v, const real5d &w,
                                         const real5d &dens,
                                         const real5d &coriolishz, int is,
                                         int js, int ks, int i, int j, int k,
                                         int n) const {
    SArray<real, 1, ndims> hv;
    // Need to subtract 1 here since d00(i,k) corresponds to p11(i,k)
    hv(0) = compute_hvxz_bottom(dens, is, js, ks, i, j, k - 1, n);
    if (ndims > 1) {
      hv(1) = compute_hvyz_bottom(dens, is, js, ks, i, j, k - 1, n);
    }

    SArray<real, 1, ndims> zetahz;
    compute_zetahz(zetahz, v, w, is, js, ks, i, j, k - 1, n);

    // compute q0 = zeta / hv and f0 = f / hv
    for (int d = 0; d < ndims; ++d) {
      qhz(d, k + ks, j + js, i + is, n) = zetahz(d) / hv(d);
      fhz(d, k + ks, j + js, i + is, n) =
          coriolishz(d, k + ks, j + js, i + is, n) / hv(d);
    }
  }

  void YAKL_INLINE compute_qxyfxy(const real5d &qxy, const real5d &fxy,
                                  const real5d &v, const real5d &w,
                                  const real5d &dens, const real5d &coriolisxy,
                                  int is, int js, int ks, int i, int j, int k,
                                  int n) const {
    real hvxy = compute_hvxy(dens, is, js, ks, i, j, k, n);

    SArray<real, 1, 1> zeta;
    compute_D1<1>(zeta, v, is, js, ks, i, j, k, n);
    qxy(0, k + ks, j + js, i + is, n) = zeta(0) / hvxy;
  }

  pvpe YAKL_INLINE compute_PVPE(const real5d &v, const real5d &w,
                                const real5d &dens, const real5d &coriolishz,
                                int is, int js, int ks, int i, int j, int k,
                                int n) const {
    pvpe vals;
    // No subtraction here since this is called on primal cells p11
    SArray<real, 1, ndims> etahz;
    compute_etahz(etahz, v, w, coriolishz, is, js, ks, i, j, k, n);
    real hv = compute_hvxz(dens, is, js, ks, i, j, k, n);
    real q0 = etahz(0) / hv;
    vals.pv = etahz(0);
    vals.pe = 0.5_fp * etahz(0) * q0;
    return vals;
  }

  pvpe YAKL_INLINE compute_PVPE_top(const real5d &v, const real5d &w,
                                    const real5d &dens,
                                    const real5d &coriolishz, int is, int js,
                                    int ks, int i, int j, int k, int n) const {
    pvpe vals;
    // No subtraction here since this is called on primal cells p11
    SArray<real, 1, ndims> etahz;
    compute_etahz(etahz, v, w, coriolishz, is, js, ks, i, j, k, n);
    real hv = compute_hvxz_top(dens, is, js, ks, i, j, k, n);
    real q0 = etahz(0) / hv;
    vals.pv = etahz(0);
    vals.pe = 0.5_fp * etahz(0) * q0;
    return vals;
  }

  pvpe YAKL_INLINE compute_PVPE_bottom(const real5d &v, const real5d &w,
                                       const real5d &dens,
                                       const real5d &coriolishz, int is, int js,
                                       int ks, int i, int j, int k,
                                       int n) const {
    pvpe vals;
    // No subtraction here since this is called on primal cells p11
    SArray<real, 1, ndims> etahz;
    compute_etahz(etahz, v, w, coriolishz, is, js, ks, i, j, k, n);
    real hv = compute_hvxz_bottom(dens, is, js, ks, i, j, k, n);
    real q0 = etahz(0) / hv;
    vals.pv = etahz(0);
    vals.pe = 0.5_fp * etahz(0) * q0;
    return vals;
  }
};
