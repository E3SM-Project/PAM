#pragma once

#include "common.h"
#include "hodge_star.h"
#include "thermo.h"
#include "variableset.h"

// ADD p-variants

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
  void YAKL_INLINE compute_dHsdx(const real3d &B, const real3d &dens,
                                 const real3d &geop, int ks, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
    compute_H2bar_vert<1, vert_diff_ord>(geop0, geop, this->primal_geometry,
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

// SHOULD BE MERGABLE INTO A SINGLE CLASS WITH INDEXING FOR RHO/RHOD?
// OR AT LEAST VARIOUS FLAGS?
// IE CE MODEL HAS CHOICES OF THERMO
// MCE MODEL HAS CHOICES OF PREDICTED VARS IE RHO VS RHOD, AND ALSO THERMO

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

  // THIS NEEDS FIXING, BUT SHOULD BE COMPLETELY GENERAL NOW!
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

    if (addmode == ADD_MODE::REPLACE) {
      B(0, k + ks, j + js, i + is, n) =
          fac * (geop0(0) + U + p * alpha - entropic_var * generalized_Exner +
                 qv * (generalized_chemical_potential_d -
                       generalized_chemical_potential_v) +
                 ql * (generalized_chemical_potential_d -
                       generalized_chemical_potential_l) +
                 qi * (generalized_chemical_potential_d -
                       generalized_chemical_potential_i));
      B(1, k + ks, j + js, i + is, n) = fac * generalized_Exner;

      // THESE INDICES ARE BROKEN!
      B(varset.dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n) =
          fac *
          (generalized_chemical_potential_v - generalized_chemical_potential_d);
      if (varset.liquid_found) {
        B(varset.dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n) =
            fac * (generalized_chemical_potential_l -
                   generalized_chemical_potential_d);
      }
      if (varset.ice_found) {
        B(varset.dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n) =
            fac * (generalized_chemical_potential_i -
                   generalized_chemical_potential_d);
      }
    } else if (addmode == ADD_MODE::ADD) {
      B(0, k + ks, j + js, i + is, n) +=
          fac * (geop0(0) + U + p * alpha - entropic_var * generalized_Exner +
                 qv * (generalized_chemical_potential_d -
                       generalized_chemical_potential_v) +
                 ql * (generalized_chemical_potential_d -
                       generalized_chemical_potential_l) +
                 qi * (generalized_chemical_potential_d -
                       generalized_chemical_potential_i));
      B(1, k + ks, j + js, i + is, n) += fac * generalized_Exner;

      // THESE INDICES ARE BROKEN!
      B(varset.dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n) +=
          fac *
          (generalized_chemical_potential_v - generalized_chemical_potential_d);
      if (varset.liquid_found) {
        B(varset.dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n) +=
            fac * (generalized_chemical_potential_l -
                   generalized_chemical_potential_d);
      }
      if (varset.ice_found) {
        B(varset.dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n) +=
            fac * (generalized_chemical_potential_i -
                   generalized_chemical_potential_d);
      }
    }
  }
};

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
