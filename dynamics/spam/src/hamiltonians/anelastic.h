#pragma once

#include "common.h"
#include "hodge_star.h"
#include "thermo.h"

// AN/MAN variants
// HOW DO WE TREAT DHDX, ETC. HERE SINCE WE HAVE A NEW ARGUMENT: P^PRIME?
// P^PRIME SOLVE IS NEEDED FOR DHDX in Hs and also for computing IE in Hs
// SHOULD THESE DETERMINE LHS AND RHS OF P EQN? YES, IT IS A CONSTRAINT BASED
// EQN THAT DEPENDS ON THERMO AND VARIABLES.... MAYBE PPRIME LIVES IN DHDX LIKE
// HS/GEOP PROPERLY SHOULD? YES ACTUALLY, THIS MAKES THE MOST SENSE...
// class Functional_PVPE_AN {};
// class Hamiltonian_Hk_AN {};
// class Hamiltonian_AN_Hs {};
// class Hamiltonian_MAN_Hs {};

#ifdef _AN
class Hamiltonian_AN_Hs {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  bool is_initialized;
  ThermoPotential thermo;
  VariableSet varset;
  real g;

  Hamiltonian_AN_Hs() { this->is_initialized = false; }

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

  void set_parameters(real gin) { this->g = gin; }

  real YAKL_INLINE compute_PE(const real5d &dens, const real5d &geop, int is,
                              int js, int ks, int i, int j, int k,
                              int n) const {
    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif
    const real rho =
        varset.get_total_density(varset.reference_state.dens.data, k, ks, n);

    return rho * geop0(0);
  }

  real YAKL_INLINE compute_IE(const real5d &dens, int is, int js, int ks, int i,
                              int j, int k, int n) const {
    const real refrho =
        1.0_fp / varset.get_alpha(varset.reference_state.dens.data, k, ks, n);
    const real refentropic_var =
        varset.get_entropic_var(varset.reference_state.dens.data, k, ks, n);
    const real refp = thermo.solve_p(refrho, refentropic_var, 0, 0, 0, 0);

    const real entropic_var =
        varset.get_entropic_var(dens, k, j, i, ks, js, is, n);

    const real H = thermo.compute_H(refp, entropic_var, 0, 0, 0, 0);

    const real rho =
        varset.get_total_density(varset.reference_state.dens.data, k, ks, n);
    return rho * H;
  }

  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx(const real5d &B, const real5d &dens,
                                 const real5d &geop, int is, int js, int ks,
                                 int i, int j, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    real refrho =
        1.0_fp / varset.get_alpha(varset.reference_state.dens.data, k, ks, n);
    real refentropic_var =
        varset.get_entropic_var(varset.reference_state.dens.data, k, ks, n);
    real refp = thermo.solve_p(refrho, refentropic_var, 0, 0, 0, 0);

    real entropic_var = varset.get_entropic_var(dens, k, j, i, ks, js, is, n);

    real H = thermo.compute_H(refp, entropic_var, 0, 0, 0, 0);
    real generalized_Exner =
        thermo.compute_dHdentropic_var(refp, entropic_var, 0, 0, 0, 0);

    if (addmode == ADD_MODE::REPLACE) {
      B(varset.active_id_mass, k + ks, j + js, i + is, n) =
          fac * (geop0(0) + H - entropic_var * generalized_Exner);
      B(varset.active_id_entr, k + ks, j + js, i + is, n) =
          fac * generalized_Exner;
    } else if (addmode == ADD_MODE::ADD) {
      B(varset.active_id_mass, k + ks, j + js, i + is, n) +=
          fac * (geop0(0) + H - entropic_var * generalized_Exner);
      B(varset.active_id_entr, k + ks, j + js, i + is, n) +=
          fac * generalized_Exner;
    }
  }

  // reference state version
  template <ADD_MODE addmode = ADD_MODE::REPLACE>
  void YAKL_INLINE compute_dHsdx(const real3d &B, const real3d &dens,
                                 const real3d &geop, int ks, int k, int n,
                                 real fac = 1._fp) const {

    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, ks, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, ks, k, n);
#endif

    real refrho =
        1.0_fp / varset.get_alpha(varset.reference_state.dens.data, k, ks, n);
    real refentropic_var =
        varset.get_entropic_var(varset.reference_state.dens.data, k, ks, n);
    real refp = thermo.solve_p(refrho, refentropic_var, 0, 0, 0, 0);

    real entropic_var = varset.get_entropic_var(dens, k, ks, n);

    real H = thermo.compute_H(refp, entropic_var, 0, 0, 0, 0);
    real generalized_Exner =
        thermo.compute_dHdentropic_var(refp, entropic_var, 0, 0, 0, 0);

    if (addmode == ADD_MODE::REPLACE) {
      B(varset.active_id_mass, k + ks, n) =
          fac * (geop0(0) + H - entropic_var * generalized_Exner);
      B(varset.active_id_entr, k + ks, n) = fac * generalized_Exner;
    } else if (addmode == ADD_MODE::ADD) {
      B(varset.active_id_mass, k + ks, n) +=
          fac * (geop0(0) + H - entropic_var * generalized_Exner);
      B(varset.active_id_entr, k + ks, n) += fac * generalized_Exner;
    }
  }
};
#endif

#ifdef _MAN
class Hamiltonian_MAN_Hs {
public:
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;
  bool is_initialized;
  ThermoPotential thermo;
  VariableSet varset;
  real g;

  using VS = VariableSet;

  Hamiltonian_MAN_Hs() { this->is_initialized = false; }

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

  void set_parameters(real gin) { this->g = gin; }

  real YAKL_INLINE compute_PE(const real5d &dens, const real5d &geop, int is,
                              int js, int ks, int i, int j, int k,
                              int n) const {
    SArray<real, 1, 1> geop0;
#ifdef _EXTRUDED
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif
    const real rho =
        varset.get_total_density(varset.reference_state.dens.data, k, ks, n);

    return rho * geop0(0);
  }

  real YAKL_INLINE compute_IE(const real5d &dens, int is, int js, int ks, int i,
                              int j, int k, int n) const {
    const real refrho =
        1.0_fp / varset.get_alpha(varset.reference_state.dens.data, k, ks, n);
    const real refentropic_var =
        varset.get_entropic_var(varset.reference_state.dens.data, k, ks, n);
    const real refqv =
        varset.get_qv(varset.reference_state.dens.data, k, ks, n);
    const real refp =
        thermo.solve_p(refrho, refentropic_var, 1 - refqv, refqv, 0, 0);

    const real entropic_var =
        varset.get_entropic_var(dens, k, j, i, ks, js, is, n);

    const real qd = varset.get_qd(dens, k, j, i, ks, js, is, n);
    const real qv = varset.get_qv(dens, k, j, i, ks, js, is, n);
    real ql = 0.0_fp;
    real qi = 0.0_fp;
    if (varset.liquid_found) {
      ql = varset.get_ql(dens, k, j, i, ks, js, is, n);
    }
    if (varset.ice_found) {
      qi = varset.get_qi(dens, k, j, i, ks, js, is, n);
    }

    const real H = thermo.compute_H(refp, entropic_var, qd, qv, ql, qi);

    const real rho =
        varset.get_total_density(varset.reference_state.dens.data, k, ks, n);
    return rho * H;
  }

  void YAKL_INLINE
  compute_dHsdx(SArray<real, 1, VS::ndensity_active> &B,
                const SArray<real, 1, VS::ndensity_dycore + VS::nmoist> &q,
                real geop) const {
    const real refp = q(0);
    const real entropic_var = q(1);
    const real qd = q(2);
    const real qv = q(3);
    const real ql = q(4);
    const real qi = q(5);

    real H = thermo.compute_H(refp, entropic_var, qd, qv, ql, qi);
    real generalized_Exner =
        thermo.compute_dHdentropic_var(refp, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_d =
        thermo.compute_dHdqd(refp, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_v =
        thermo.compute_dHdqv(refp, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_l =
        thermo.compute_dHdql(refp, entropic_var, qd, qv, ql, qi);
    real generalized_chemical_potential_i =
        thermo.compute_dHdqi(refp, entropic_var, qd, qv, ql, qi);

    B(varset.active_id_mass) = geop + H - entropic_var * generalized_Exner +
                               qv * (generalized_chemical_potential_d -
                                     generalized_chemical_potential_v) +
                               ql * (generalized_chemical_potential_d -
                                     generalized_chemical_potential_l) +
                               qi * (generalized_chemical_potential_d -
                                     generalized_chemical_potential_i);
    B(varset.active_id_entr) = generalized_Exner;

    if (!ThermoPotential::moist_species_decouple_from_dynamics) {
      B(varset.active_id_vap) =
          generalized_chemical_potential_v - generalized_chemical_potential_d;
      if (varset.liquid_found) {
        B(varset.active_id_liq) =
            generalized_chemical_potential_l - generalized_chemical_potential_d;
      }
      if (varset.ice_found) {
        B(varset.active_id_ice) =
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
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, is, js, ks, i,
        j, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, is, js, ks, i, j, k, n);
#endif

    const real refrho =
        1.0_fp / varset.get_alpha(varset.reference_state.dens.data, k, ks, n);
    const real refentropic_var =
        varset.get_entropic_var(varset.reference_state.dens.data, k, ks, n);
    const real refqv =
        varset.get_qv(varset.reference_state.dens.data, k, ks, n);
    const real refp =
        thermo.solve_p(refrho, refentropic_var, 1 - refqv, refqv, 0, 0);

    SArray<real, 1, VS::ndensity_active> l_B;
    SArray<real, 1, VS::ndensity_dycore + VS::nmoist> l_q;

    l_q(0) = refp;
    l_q(1) = varset.get_entropic_var(dens, k, j, i, ks, js, is, n);
    l_q(2) = varset.get_qd(dens, k, j, i, ks, js, is, n);
    l_q(3) = varset.get_qv(dens, k, j, i, ks, js, is, n);
    l_q(4) = varset.liquid_found ? varset.get_ql(dens, k, j, i, ks, js, is, n)
                                 : 0.0_fp;
    l_q(5) =
        varset.ice_found ? varset.get_qi(dens, k, j, i, ks, js, is, n) : 0.0_fp;

    compute_dHsdx(l_B, l_q, geop0(0));

    for (int d = 0; d < VS::ndensity_active; ++d) {
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
#ifdef _EXTRUDED
    compute_Hn1bar<1, diff_ord, vert_diff_ord>(
        geop0, geop, this->primal_geometry, this->dual_geometry, ks, k, n);
#else
    compute_H2bar<1, diff_ord>(geop0, geop, this->primal_geometry,
                               this->dual_geometry, ks, k, n);
#endif

    const real refrho =
        1.0_fp / varset.get_alpha(varset.reference_state.dens.data, k, ks, n);
    const real refentropic_var =
        varset.get_entropic_var(varset.reference_state.dens.data, k, ks, n);
    const real refqv =
        varset.get_qv(varset.reference_state.dens.data, k, ks, n);
    const real refp =
        thermo.solve_p(refrho, refentropic_var, 1 - refqv, refqv, 0, 0);

    SArray<real, 1, VS::ndensity_active> l_B;
    SArray<real, 1, VS::ndensity_dycore + VS::nmoist> l_q;

    l_q(0) = refp;
    l_q(1) = varset.get_entropic_var(dens, k, ks, n);
    l_q(2) = varset.get_qd(dens, k, ks, n);
    l_q(3) = varset.get_qv(dens, k, ks, n);
    l_q(4) = 0.0_fp;
    l_q(5) = 0.0_fp;

    compute_dHsdx(l_B, l_q, geop0(0));

    for (int d = 0; d < VS::ndensity_active; ++d) {
      if (addmode == ADD_MODE::REPLACE) {
        B(d, k + ks, n) = fac * l_B(d);
      } else if (addmode == ADD_MODE::ADD) {
        B(d, k + ks, n) += fac * l_B(d);
      }
    }
  }
};
#endif
