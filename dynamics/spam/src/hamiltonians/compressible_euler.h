#pragma once

#include "common.h"
#include "hodge_star.h"
#include "thermo.h"
#include "variableset.h"

real YAKL_INLINE gamma_avg(real a, real b, real gamma) {
  const real f = (a - b) / (a + b);
  const real v = f * f;
  if (v < 1e-4_fp) {
    const real c1 = 1._fp / 6._fp * (gamma - 1) * (gamma - 2);
    const real c2 = 1._fp / 20._fp * (gamma - 3) * (gamma - 4);
    const real c3 = 1._fp / 42._fp * (gamma - 5) * (gamma - 6);
    const real x = std::pow(0.5_fp * (a + b), gamma - 1);
    return x * (1 + c1 * v * (1 + c2 * v * (1 + c3 * v)));
  } else {
    return (std::pow(a, gamma) - std::pow(b, gamma)) / (gamma * (a - b));
  }
}

// ADD p-variants

#ifdef _CE
class Hamiltonian_CE_Hs {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  bool is_initialized;
  ThermoPotential thermo;
  VariableSet varset;

  Hamiltonian_CE_Hs() { this->is_initialized = false; }

  void initialize(const ThermoPotential &thermodynamics,
                  const VariableSet &variableset,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->thermo = thermodynamics;
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->is_initialized = true;
    this->varset = variableset;
  }

  void set_parameters(real gin) {}

  real YAKL_INLINE compute_PE(const real5d &dens, const real5d &geop, int is,
                              int js, int ks, int i, int j, int k,
                              int n) const {
    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    return dens(0, k + ks, j + js, i + is, n) * geop0(0);
  }

  real YAKL_INLINE compute_IE(const real5d &dens, int is, int js, int ks, int i,
                              int j, int k, int n) const {
    // SArray<real,1,2> dens0;
    // #ifdef _EXTRUDED
    // compute_Iext<2, diff_ord, vert_diff_ord>(dens0, dens,
    // *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
    // #else
    // compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry,
    // *this->dual_geometry, is, js, ks, i, j, k, n); #endif

    real alpha = varset.get_alpha(dens, k, j, i, ks, js, is, n);
    real entropic_var = varset.get_entropic_var(dens, k, j, i, ks, js, is, n);
    return dens(0, k + ks, j + js, i + is, n) *
           thermo.compute_U(alpha, entropic_var, 0._fp, 0._fp, 0._fp, 0._fp);
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx(const real5d &B, const real5d &dens,
                                 const real5d &geop, int is, int js, int ks,
                                 int i, int j, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    // SArray<real,1,2> dens0;
    // #ifdef _EXTRUDED
    // compute_Iext<2, diff_ord, vert_diff_ord>(dens0, dens,
    // *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
    // #else
    // compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry,
    // *this->dual_geometry, is, js, ks, i, j, k, n); #endif

    // real alpha = 1.0_fp / dens0(0);
    // real entropic_var = dens(1)/dens(0);

    // real alpha = 1.0_fp / dens(0,k+ks,j+js,i+is,n);
    // real entropic_var = dens(1,k+ks,j+js,i+is,n)/dens(0,k+ks,j+js,i+is,n);

    real alpha = varset.get_alpha(dens, k, j, i, ks, js, is, n);
    real entropic_var = varset.get_entropic_var(dens, k, j, i, ks, js, is, n);

    real U = thermo.compute_U(alpha, entropic_var, 0, 0, 0, 0);
    real p = -thermo.compute_dUdalpha(alpha, entropic_var, 0, 0, 0, 0);
    real generalized_Exner =
        thermo.compute_dUdentropic_var(alpha, entropic_var, 0, 0, 0, 0);

    if (addmode == ADD_MODE::REPLACE) {
      B(0, k + ks, j + js, i + is, n) =
          fac * (geop0(0) + U + p * alpha - entropic_var * generalized_Exner);
      B(1, k + ks, j + js, i + is, n) = fac * generalized_Exner;
    } else if (addmode == ADD_MODE::ADD) {
      B(0, k + ks, j + js, i + is, n) +=
          fac * (geop0(0) + U + p * alpha - entropic_var * generalized_Exner);
      B(1, k + ks, j + js, i + is, n) += fac * generalized_Exner;
    }
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx_two_point(const real5d &B, const real5d &dens1,
                                           const real5d &dens2,
                                           const real5d &geop, int is, int js,
                                           int ks, int i, int j, int k, int n,
                                           real fac = 1._fp) const {
    // hacky way to prevent compilation errors for unimplemented thermo variants
    if constexpr (!si_compute_functional_derivatives_quadrature) {
      // dispatch based on thermo
      compute_dHsdx_two_point<addmode>(thermo, B, dens1, dens2, geop, is, js,
                                       ks, i, j, k, n, fac);
    }
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx_two_point(IdealGas_Pottemp, const real5d &B,
                                           const real5d &dens1,
                                           const real5d &dens2,
                                           const real5d &geop, int is, int js,
                                           int ks, int i, int j, int k, int n,
                                           real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    real Cpd = thermo.cst.Cpd;
    real Cvd = thermo.cst.Cvd;
    real Rd = thermo.cst.Rd;
    real pr = thermo.cst.pr;
    real gamma_d = thermo.cst.gamma_d;

    real Tht1 = dens1(1, k + ks, j + js, i + is, n) /
                dual_geometry.get_area_11entity(k + ks, j + js, i + is);
    real Tht2 = dens2(1, k + ks, j + js, i + is, n) /
                dual_geometry.get_area_11entity(k + ks, j + js, i + is);

    real generalized_Exner =
        Cpd * std::pow(Rd / pr, gamma_d - 1) * gamma_avg(Tht1, Tht2, gamma_d);

    if (addmode == ADD_MODE::REPLACE) {
      B(0, k + ks, j + js, i + is, n) = fac * geop0(0);
      B(1, k + ks, j + js, i + is, n) = fac * generalized_Exner;
    } else if (addmode == ADD_MODE::ADD) {
      B(0, k + ks, j + js, i + is, n) += fac * geop0(0);
      B(1, k + ks, j + js, i + is, n) += fac * generalized_Exner;
    }
  }

  // reference state version
  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx(const real3d &B, const real3d &dens,
                                 const real3d &geop, int ks, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
    compute_H2bar_ext<1, vert_diff_ord>(geop0, geop, this->primal_geometry,
                                        this->dual_geometry, ks, k, n);

    real alpha = varset.get_alpha(dens, k, ks, n);
    real entropic_var = varset.get_entropic_var(dens, k, ks, n);

    real U = thermo.compute_U(alpha, entropic_var, 0, 0, 0, 0);
    real p = -thermo.compute_dUdalpha(alpha, entropic_var, 0, 0, 0, 0);
    real generalized_Exner =
        thermo.compute_dUdentropic_var(alpha, entropic_var, 0, 0, 0, 0);

    if (addmode == ADD_MODE::REPLACE) {
      B(0, k, n) =
          fac * (geop0(0) + U + p * alpha - entropic_var * generalized_Exner);
      B(1, k, n) = fac * generalized_Exner;
    } else if (addmode == ADD_MODE::ADD) {
      B(0, k, n) +=
          fac * (geop0(0) + U + p * alpha - entropic_var * generalized_Exner);
      B(1, k, n) += fac * generalized_Exner;
    }
  }
};
#endif

// SHOULD BE MERGABLE INTO A SINGLE CLASS WITH INDEXING FOR RHO/RHOD?
// OR AT LEAST VARIOUS FLAGS?
// IE CE MODEL HAS CHOICES OF THERMO
// MCE MODEL HAS CHOICES OF PREDICTED VARS IE RHO VS RHOD, AND ALSO THERMO

#ifdef _MCErho
class Hamiltonian_MCE_Hs {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  bool is_initialized;
  ThermoPotential thermo;
  VariableSet varset;

  Hamiltonian_MCE_Hs() { this->is_initialized = false; }

  void initialize(const ThermoPotential &thermodynamics,
                  const VariableSet &variableset,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom) {
    this->thermo = thermodynamics;
    this->primal_geometry = primal_geom;
    this->dual_geometry = dual_geom;
    this->is_initialized = true;
    this->varset = variableset;
  }

  void set_parameters(real gin) {}

  real YAKL_INLINE compute_PE(const real5d &dens, const real5d &geop, int is,
                              int js, int ks, int i, int j, int k,
                              int n) const {
    SArray<real, 1, 1> geop0;

#ifdef _EXTRUDED
    compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    return dens(0, k + ks, j + js, i + is, n) * geop0(0);
  }

  real YAKL_INLINE compute_IE(const real5d &dens, int is, int js, int ks, int i,
                              int j, int k, int n) const {

    // SArray<real,1,5> dens0;
    // #ifdef _EXTRUDED
    // compute_Iext<5, diff_ord, vert_diff_ord>(dens0, dens,
    // *this->primal_geometry, *this->dual_geometry, is, js, ks, i, j, k, n);
    // #else
    // compute_I<5, diff_ord>(dens0, dens, *this->primal_geometry,
    // *this->dual_geometry, is, js, ks, i, j, k, n); #endif

    real alpha = varset.get_alpha(dens, k, j, i, ks, js, is, n);
    real entropic_var = varset.get_entropic_var(dens, k, j, i, ks, js, is, n);
    real qd = varset.get_qd(dens, k, j, i, ks, js, is, n);
    real qv = varset.get_qv(dens, k, j, i, ks, js, is, n);
    real ql = 0.0_fp;
    real qi = 0.0_fp;
    if (varset.liquid_found) {
      ql = varset.get_ql(dens, k, j, i, ks, js, is, n);
    }
    if (varset.ice_found) {
      qi = varset.get_qi(dens, k, j, i, ks, js, is, n);
    }
    return varset.get_total_density(dens, k, j, i, ks, js, is, n) *
           thermo.compute_U(alpha, entropic_var, qd, qv, ql, qi);
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx_two_point(const real5d &B, const real5d &dens1,
                                           const real5d &dens2,
                                           const real5d &geop, int is, int js,
                                           int ks, int i, int j, int k, int n,
                                           real fac = 1._fp) const {
    // hacky way to prevent compilation errors for unimplemented thermo variants
    if constexpr (!si_compute_functional_derivatives_quadrature) {
      // dispatch based on thermo
      compute_dHsdx_two_point<addmode>(thermo, B, dens1, dens2, geop, is, js,
                                       ks, i, j, k, n, fac);
    }
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx_two_point(ConstantKappa_VirtualPottemp,
                                           const real5d &B, const real5d &dens1,
                                           const real5d &dens2,
                                           const real5d &geop, int is, int js,
                                           int ks, int i, int j, int k, int n,
                                           real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    real Cpd = thermo.cst.Cpd;
    real Cvd = thermo.cst.Cvd;
    real Rd = thermo.cst.Rd;
    real pr = thermo.cst.pr;
    real gamma_d = thermo.cst.gamma_d;

    real Tht1 = dens1(1, k + ks, j + js, i + is, n) /
                dual_geometry.get_area_11entity(k + ks, j + js, i + is);
    real Tht2 = dens2(1, k + ks, j + js, i + is, n) /
                dual_geometry.get_area_11entity(k + ks, j + js, i + is);

    real generalized_Exner =
        Cpd * std::pow(Rd / pr, gamma_d - 1) * gamma_avg(Tht1, Tht2, gamma_d);

    if (addmode == ADD_MODE::REPLACE) {
      B(0, k + ks, j + js, i + is, n) = fac * geop0(0);
      B(1, k + ks, j + js, i + is, n) = fac * generalized_Exner;
    } else if (addmode == ADD_MODE::ADD) {
      B(0, k + ks, j + js, i + is, n) += fac * geop0(0);
      B(1, k + ks, j + js, i + is, n) += fac * generalized_Exner;
    }

    // assumes that tracers_decouple_from_dynamics == true !
    // TODO: how to statically assert this ?
  }

  void YAKL_INLINE compute_dHsdx(
      SArray<real, 1, ndensity_B> &B,
      const SArray<real, 1, ndensity_dycore + nmoist> &q, real geop) const {
    const real alpha = q(0);
    const real entropic_var = q(1);
    const real qd = q(2);
    const real qv = q(3);
    const real ql = q(4);
    const real qi = q(5);

    real U = thermo.compute_U(alpha, entropic_var, qd, qv, ql, qi);
    real p = -thermo.compute_dUdalpha(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_Exner =
        thermo.compute_dUdentropic_var(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_d =
        thermo.compute_dUdqd(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_v =
        thermo.compute_dUdqv(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_l =
        thermo.compute_dUdql(alpha, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_i =
        thermo.compute_dUdqi(alpha, entropic_var, qd, qv, ql, qi);

    B(0) = geop + U + p * alpha - entropic_var * generalized_Exner +
           qv * (generalized_chemical_potential_d -
                 generalized_chemical_potential_v) +
           ql * (generalized_chemical_potential_d -
                 generalized_chemical_potential_l) +
           qi * (generalized_chemical_potential_d -
                 generalized_chemical_potential_i);
    B(1) = generalized_Exner;

    if (!tracers_decouple_from_dynamics) {
      B(varset.dm_id_vap + ndensity_nophysics) =
          generalized_chemical_potential_v - generalized_chemical_potential_d;
      if (varset.liquid_found) {
        B(varset.dm_id_liq + ndensity_nophysics) =
            generalized_chemical_potential_l - generalized_chemical_potential_d;
      }
      if (varset.ice_found) {
        B(varset.dm_id_ice + ndensity_nophysics) =
            generalized_chemical_potential_i - generalized_chemical_potential_d;
      }
    }
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx(const real5d &B, const real5d &dens,
                                 const real5d &geop, int is, int js, int ks,
                                 int i, int j, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;

#ifdef _EXTRUDED
    compute_H2bar_ext<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    SArray<real, 1, ndensity_B> l_B;
    SArray<real, 1, ndensity_dycore + nmoist> l_q;

    l_q(0) = varset.get_alpha(dens, k, j, i, ks, js, is, n);
    l_q(1) = varset.get_entropic_var(dens, k, j, i, ks, js, is, n);
    l_q(2) = varset.get_qd(dens, k, j, i, ks, js, is, n);
    l_q(3) = varset.get_qv(dens, k, j, i, ks, js, is, n);
    l_q(4) = varset.liquid_found ? varset.get_ql(dens, k, j, i, ks, js, is, n)
                                 : 0.0_fp;
    l_q(5) = varset.liquid_found ? varset.get_qi(dens, k, j, i, ks, js, is, n)
                                 : 0.0_fp;

    compute_dHsdx(l_B, l_q, geop0(0));

    for (int d = 0; d < ndensity_B; ++d) {
      if (addmode == ADD_MODE::REPLACE) {
        B(d, k + ks, j + js, i + is, n) = fac * l_B(d);
      } else if (addmode == ADD_MODE::ADD) {
        B(d, k + ks, j + js, i + is, n) += fac * l_B(d);
      }
    }
  }

  // reference state version
  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx(const real3d &B, const real3d &dens,
                                 const real3d &geop, int ks, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
    compute_H2bar_ext<1, vert_diff_ord>(geop0, geop, this->primal_geometry,
                                        this->dual_geometry, ks, k, n);

    SArray<real, 1, ndensity_B> l_B;
    SArray<real, 1, ndensity_dycore + nmoist> l_q;

    l_q(0) = varset.get_alpha(dens, k, ks, n);
    l_q(1) = varset.get_entropic_var(dens, k, ks, n);
    l_q(2) = varset.get_qd(dens, k, ks, n);
    l_q(3) = varset.get_qv(dens, k, ks, n);
    l_q(4) = 0.0_fp;
    l_q(5) = 0.0_fp;

    compute_dHsdx(l_B, l_q, geop0(0));

    for (int d = 0; d < ndensity_B; ++d) {
      if (addmode == ADD_MODE::REPLACE) {
        B(d, k, n) = fac * l_B(d);
      } else if (addmode == ADD_MODE::ADD) {
        B(d, k, n) += fac * l_B(d);
      }
    }
  }
};
#endif

// //ALL BROKEN, BUT REALLY CAN MOSTLY BE ELIMINATED?
//
// class Hamiltonian_MCE_rhod_Hs {
//  public:
//    Geometry<1,1,1> *primal_geometry;
//    Geometry<1,1,1> *dual_geometry;
//    bool is_initialized;
//    ThermoPotential *thermo;
//
//     Hamiltonian_MCE_rhod_Hs() {, int n
//       this->is_initialized = false;
//  }
//
//  void initialize(ModelParameters &params, ThermoPotential &thermodynamics,
//  Geometry<1,1,1> &primal_geom, Geometry<1,1,1> &dual_geom)
//  {
//    this->thermo = &thermodynamics;
//    this->primal_geometry = &primal_geom;
//    this->dual_geometry = &dual_geom;
//    this->is_initialized = true;
//  }
//
//
//  real YAKL_INLINE compute_PE(const real5d& dens, const real5d& densfct, const
//  const real5d& geop, int is, int js, int ks, int i, int j, int k)
//  {
//
//    SArray<real,1,1> geop0;
//    compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry,
//    *this->dual_geometry, is, js, ks, i, j, k);
//
//    return (dens(0, k+ks, j+js, i+is) + densfct(0, k+ks, j+js, i+is) +
//    densfct(1, k+ks, j+js, i+is) + densfct(2, k+ks, j+js, i+is) ) * geop0(0);
//
//  }
//
//  //   alpha = 1./rho;
//  //   compute_entropic_var = entropic_density/rho;
//  //   qd = rho_d/rho;
//  //   qv = rho_v/rho;
//  //   ql = rho_l/rho;
//  //   qi = rho_i/rho;
//  // with rho = rho_d + rho_v + rho_l + rho_i
//      real YAKL_INLINE compute_alpha(const real5d& dens0, const real5d&
//      densfct0, int is, int js, int ks, int i, int j, int k, int n) {
//        return 1./(dens0(0, k+ks, j+js, i+is) + densfct0(0, k+ks, j+js, i+is)
//        + densfct0(1, k+ks, j+js, i+is) + densfct0(2, k+ks, j+js, i+is));
//      }
//
//      real YAKL_INLINE compute_entropic_var(const real5d& dens0, const real5d&
//      densfct0, int is, int js, int ks, int i, int j, int k, int n) {
//        return dens0(1, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) +
//        densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) +
//        densfct0(2, k+ks, j+js, i+is));
//      }
//      real YAKL_INLINE compute_qd(const real5d& dens0, const real5d& densfct0,
//      int is, int js, int ks, int i, int j, int k, int n) {
//        return (dens0(0, k+ks, j+js, i+is))/(dens0(0, k+ks, j+js, i+is) +
//        densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) +
//        densfct0(2, k+ks, j+js, i+is));
//      }
//      real YAKL_INLINE compute_qv(const real5d& dens0, const real5d& densfct0,
//      int is, int js, int ks, int i, int j, int k, int n) {
//        return densfct0(0, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) +
//        densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) +
//        densfct0(2, k+ks, j+js, i+is));
//      }
//      real YAKL_INLINE compute_ql(const real5d& dens0, const real5d& densfct0,
//      int is, int js, int ks, int i, int j, int k, int n) {
//        return densfct0(1, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) +
//        densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) +
//        densfct0(2, k+ks, j+js, i+is));
//      }
//      real YAKL_INLINE compute_qi(const real5d& dens0, const real5d& densfct0,
//      int is, int js, int ks, int i, int j, int k, int n) {
//        return densfct0(2, k+ks, j+js, i+is)/(dens0(0, k+ks, j+js, i+is) +
//        densfct0(0, k+ks, j+js, i+is) + densfct0(1, k+ks, j+js, i+is) +
//        densfct0(2, k+ks, j+js, i+is));
//      }
//
//
//  real YAKL_INLINE compute_IE(const real5d& dens, const real5d& densfct, int
//  is, int js, int ks, int i, int j, int k, int n)
//  {
//
//    SArray<real,1,2> dens0;
//    SArray<real,1,3> densfct0;
//    compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry,
//    *this->dual_geometry, is, js, ks, i, j, k); compute_I<3,
//    diff_ord>(densfct0, densfct, *this->primal_geometry, *this->dual_geometry,
//    is, js, ks, i, j, k); real rho = dens0(0) + densfct(0) + densfct(1) +
//    densfct(2); real alpha = 1./rho; real entropic_var = dens0(1)/rho; real qd
//    = dens0(0) / rho; real qv = densfct(0) / rho; real ql = densfct(1) / rho;
//    real qi = densfct(2) / rho;
//    return rho * thermo.compute_U(alpha, entropic_var, qd, qv, ql, qi);
//  }
//
//   void YAKL_INLINE compute_dHsdx(const real5d& B, const real5d& Bfct, const
//   real5d& dens0, const real5d& densfct0, const real5d& geop, int is, int js,
//   int ks, int i, int j, int k, int n)
//   {
//
//     SArray<real,1,1> geop0;
//     compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry,
//     *this->dual_geometry, is, js, ks, i, j, k);
//
//     real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
//     real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i,
//     j, k); real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k); real
//     qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k); real ql =
//     compute_ql(dens0, densfct0, is, js, ks, i, j, k); real qi =
//     compute_qi(dens0, densfct0, is, js, ks, i, j, k);
//
//     real U = thermo.compute_U(alpha, entropic_var, qd, qv, ql, qi);
//     real p = -thermo.compute_dUdalpha(alpha, entropic_var, qd, qv, ql, qi);
//     real generalized_Exner = thermo.compute_dUdentropic_var(alpha,
//     entropic_var, qd, qv, ql, qi); real generalized_chemical_potential_d =
//     thermo.compute_dUdqd(alpha, entropic_var, qd, qv, ql, qi); real
//     generalized_chemical_potential_v = thermo.compute_dUdqv(alpha,
//     entropic_var, qd, qv, ql, qi); real generalized_chemical_potential_l =
//     thermo.compute_dUdql(alpha, entropic_var, qd, qv, ql, qi); real
//     generalized_chemical_potential_i = thermo.compute_dUdqi(alpha,
//     entropic_var, qd, qv, ql, qi);
//
//     real QNterm = qd * generalized_chemical_potential_d + qv *
//     generalized_chemical_potential_v + qi * generalized_chemical_potential_i
//     + qi * generalized_chemical_potential_i; B(0, k+ks, j+js, i+is) =
//     geop0(0) + U + p * alpha - entropic_var * generalized_Exner -
//     generalized_chemical_potential_d - QNterm; Bfct(0, k+ks, j+js, i+is) =
//     geop0(0) + U + p * alpha - entropic_var * generalized_Exner -
//     generalized_chemical_potential_v - QNterm; Bfct(1, k+ks, j+js, i+is) =
//     geop0(0) + U + p * alpha - entropic_var * generalized_Exner -
//     generalized_chemical_potential_l - QNterm; Bfct(2, k+ks, j+js, i+is) =
//     geop0(0) + U + p * alpha - entropic_var * generalized_Exner -
//     generalized_chemical_potential_i - QNterm; B(1, k+ks, j+js, i+is) =
//     generalized_Exner;
//   }
// };

// p-variants
// class Hamiltonian_CE_p_Hs : public Hamiltonian_CE_Hs
// {
//   real YAKL_INLINE compute_IE(const real5d& dens, const real5d& densfct, int
//   is, int js, int ks, int i, int j, int k)
//   {
//     SArray<real,2> dens0;
//     compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry,
//     *this->dual_geometry, is, js, ks, i, j, k); real entropic_var =
//     dens0(1)/dens0(0); real p = thermo.solve_p(dens0(0), entropic_var, 0, 0,
//     0, 0); return dens(0, k+ks, j+js, i+is) * thermo.compute_H(p,
//     entropic_var, 0, 0, 0, 0) - p;
//   }
//
//    void YAKL_INLINE compute_dHsdx(const real5d& B, const real5d& Bfct, const
//    real5d& dens0, const real5d& densfct0, const real5d& geop, int is, int js,
//    int ks, int i, int j, int k)
//    {
//
//      SArray<real,1> geop0;
//      compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry,
//      *this->dual_geometry, is, js, ks, i, j, k);
//
//      real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i,
//      j, k); real p = thermo.solve_p(dens0(0, k+ks, j+js, i+is),
//      entropic_var, 0, 0, 0, 0);
//
//      real H = thermo.compute_H(p, entropic_var, 0, 0, 0, 0);
//      real generalized_Exner = thermo.compute_dHdentropic_var(p,
//      entropic_var, 0, 0, 0, 0);
//
//      B(0, k+ks, j+js, i+is) = geop0(0) + H - entropic_var *
//      generalized_Exner; B(1, k+ks, j+js, i+is) = generalized_Exner;
//
//    }
// };
//
// class Hamiltonian_MCE_rho_p_Hs : public Hamiltonian_MCE_rho_Hs
// {
//   real YAKL_INLINE compute_IE(const real5d& dens, const real5d& densfct, int
//   is, int js, int ks, int i, int j, int k)
//   {
//
//     SArray<real,2> dens0;
//     SArray<real,3> densfct0;
//     compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry,
//     *this->dual_geometry, is, js, ks, i, j, k); compute_I<3,
//     diff_ord>(densfct0, densfct, *this->primal_geometry,
//     *this->dual_geometry, is, js, ks, i, j, k);
//
//     real entropic_var = dens0(1)/dens0(0);
//     real qd = (dens0(0) - densfct(0) - densfct(1) - densfct(2))/dens0(0);
//     real qv = densfct(0) / dens0(0);
//     real ql = densfct(1) / dens0(0);
//     real qi = densfct(2) / dens0(0);
//     real p = thermo.solve_p(dens0(0), entropic_var, qd, qv, ql, qi);
//     return dens(0, k+ks, j+js, i+is) * thermo.compute_H(p, entropic_var, qd,
//     qv, ql, qi) - p;
//   }
//
//    void YAKL_INLINE compute_dHsdx(const real5d& B, const real5d& Bfct, const
//    real5d& dens0, const real5d& densfct0, const real5d& geop, int is, int js,
//    int ks, int i, int j, int k)
//    {
//
//      SArray<real,1> geop0;
//      compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry,
//      *this->dual_geometry, is, js, ks, i, j, k);
//
//      real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i,
//      j, k); real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k); real
//      qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k); real ql =
//      compute_ql(dens0, densfct0, is, js, ks, i, j, k); real qi =
//      compute_qi(dens0, densfct0, is, js, ks, i, j, k);
//
//      real p = thermo.solve_p(dens0(0, k+ks, j+js, i+is), entropic_var, qd,
//      qv, ql, qi); real H = thermo.compute_H(p, entropic_var, qd, qv, ql,
//      qi); real generalized_Exner = thermo.compute_dHdentropic_var(p,
//      entropic_var, qd, qv, ql, qi); real generalized_chemical_potential_d =
//      thermo.compute_dHdqd(p, entropic_var, qd, qv, ql, qi); real
//      generalized_chemical_potential_v = thermo.compute_dHdqv(p,
//      entropic_var, qd, qv, ql, qi); real generalized_chemical_potential_l =
//      thermo.compute_dHdql(p, entropic_var, qd, qv, ql, qi); real
//      generalized_chemical_potential_i = thermo.compute_dHdqi(p,
//      entropic_var, qd, qv, ql, qi);
//
//      B(0, k+ks, j+js, i+is) = geop0(0) + H - entropic_var * generalized_Exner
//      + qv * (generalized_chemical_potential_d -
//      generalized_chemical_potential_v) + ql *
//      (generalized_chemical_potential_d - generalized_chemical_potential_l) +
//      qi * (generalized_chemical_potential_d -
//      generalized_chemical_potential_i); B(1, k+ks, j+js, i+is) =
//      generalized_Exner; Bfct(0, k+ks, j+js, i+is) =
//      generalized_chemical_potential_v - generalized_chemical_potential_d;
//      Bfct(1, k+ks, j+js, i+is) = generalized_chemical_potential_l -
//      generalized_chemical_potential_d; Bfct(2, k+ks, j+js, i+is) =
//      generalized_chemical_potential_i - generalized_chemical_potential_d;
//
//    }
// };
//
// class Hamiltonian_MCE_rhod_p_Hs : public Hamiltonian_MCE_rhod_Hs
// {
//   real YAKL_INLINE compute_IE(const real5d& dens, const real5d& densfct, int
//   is, int js, int ks, int i, int j, int k)
//   {
//
//     SArray<real,2> dens0;
//     SArray<real,3> densfct0;
//     compute_I<2, diff_ord>(dens0, dens, *this->primal_geometry,
//     *this->dual_geometry, is, js, ks, i, j, k); compute_I<3,
//     diff_ord>(densfct0, densfct, *this->primal_geometry,
//     *this->dual_geometry, is, js, ks, i, j, k);
//
//     real rho = dens0(0) + densfct(0) + densfct(1) + densfct(2);
//     real entropic_var = dens0(1)/rho;
//     real qd = dens0(0) / rho;
//     real qv = densfct(0) / rho;
//     real ql = densfct(1) / rho;
//     real qi = densfct(2) / rho;
//     real p = thermo.solve_p(rho, entropic_var, qd, qv, ql, qi);
//     return rho * thermo.compute_H(p, entropic_var, qd, qv, ql, qi) - p;
//   }
//
//    void YAKL_INLINE compute_dHsdx(const real5d& B, const real5d& Bfct, const
//    real5d& dens0, const real5d& densfct0, const real5d& geop, int is, int js,
//    int ks, int i, int j, int k)
//    {
//
//      SArray<real,1> geop0;
//      compute_I<1,diff_ord> (geop0, geop, *this->primal_geometry,
//      *this->dual_geometry, is, js, ks, i, j, k);
//
//      real alpha = compute_alpha(dens0, densfct0, is, js, ks, i, j, k);
//      real entropic_var = compute_entropic_var(dens0, densfct0, is, js, ks, i,
//      j, k); real qd = compute_qd(dens0, densfct0, is, js, ks, i, j, k); real
//      qv = compute_qv(dens0, densfct0, is, js, ks, i, j, k); real ql =
//      compute_ql(dens0, densfct0, is, js, ks, i, j, k); real qi =
//      compute_qi(dens0, densfct0, is, js, ks, i, j, k); real p =
//      thermo.solve_p(1./alpha, entropic_var, qd, qv, ql, qi);
//
//      real H = thermo.compute_H(p, entropic_var, qd, qv, ql, qi);
//      real generalized_Exner = thermo.compute_dHdentropic_var(p,
//      entropic_var, qd, qv, ql, qi); real generalized_chemical_potential_d =
//      thermo.compute_dHdqd(p, entropic_var, qd, qv, ql, qi); real
//      generalized_chemical_potential_v = thermo.compute_dHdqv(p,
//      entropic_var, qd, qv, ql, qi); real generalized_chemical_potential_l =
//      thermo.compute_dHdql(p, entropic_var, qd, qv, ql, qi); real
//      generalized_chemical_potential_i = thermo.compute_dHdqi(p,
//      entropic_var, qd, qv, ql, qi);
//
//      real QNterm = qd * generalized_chemical_potential_d + qv *
//      generalized_chemical_potential_v + qi * generalized_chemical_potential_i
//      + qi * generalized_chemical_potential_i; B(0, k+ks, j+js, i+is) =
//      geop0(0) + H - entropic_var * generalized_Exner -
//      generalized_chemical_potential_d - QNterm; Bfct(0, k+ks, j+js, i+is) =
//      geop0(0) + H - entropic_var * generalized_Exner -
//      generalized_chemical_potential_v - QNterm; Bfct(1, k+ks, j+js, i+is) =
//      geop0(0) + H - entropic_var * generalized_Exner -
//      generalized_chemical_potential_l - QNterm; Bfct(2, k+ks, j+js, i+is) =
//      geop0(0) + H - entropic_var * generalized_Exner -
//      generalized_chemical_potential_i - QNterm; B(1, k+ks, j+js, i+is) =
//      generalized_Exner;
//    }
// };
//
//
