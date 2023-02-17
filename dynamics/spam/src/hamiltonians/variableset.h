#pragma once

#include "DataManager.h"
#include "MultipleFields.h"
#include "common.h"
#include "geometry.h"
#include "pam_coupler.h" //Has DataManager and pam_const
#include "refstate.h"
#include "thermo.h"
using pam::PamCoupler;

// solve a system to exactly invert the velocity averaging done
// during conversion to coupler state when coupling winds
constexpr bool couple_wind_exact_inverse = false;

struct VS_SWE {
  static constexpr bool couple = false;
};
struct VS_TSWE {
  static constexpr bool couple = false;
};
struct VS_CE {
  static constexpr bool couple = false;
};
struct VS_AN {
  static constexpr bool couple = false;
};
struct VS_MAN {
  static constexpr bool couple = true;
};
struct VS_MCE_rho {
  static constexpr bool couple = true;
};
struct VS_MCE_rhod {
  static constexpr bool couple = true;
};
struct VS_CE_p {
  static constexpr bool couple = false;
};
struct VS_MCE_rhop {
  static constexpr bool couple = false;
};
struct VS_MCE_rhodp {
  static constexpr bool couple = false;
};

template <class T> class VariableSetBase {
public:
  std::string dens_name[ndensity]; // Name of each density
  std::string dens_desc[ndensity]; // Description of each density
  // bool1d                   dens_pos;       // Whether each density is
  // positive-definite
  SArray<bool, 1, ndensity>
      dens_pos; // Whether each density is positive-definite
  bool couple_wind;

  int dm_id_vap = -1;
  int dm_id_liq = -1;
  int dm_id_ice = -1;

  bool ice_found = false;
  bool liquid_found = false;

  ThermoPotential thermo;
  ReferenceState reference_state;
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  void initialize(PamCoupler &coupler, ModelParameters &params,
                  const ThermoPotential &thermodynamics,
                  const ReferenceState &refstate,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom);
  static void initialize(VariableSetBase &varset, PamCoupler &coupler,
                         ModelParameters &params,
                         const ThermoPotential &thermodynamics,
                         const ReferenceState &refstate,
                         const Geometry<Straight> &primal_geom,
                         const Geometry<Twisted> &dual_geom) {

    varset.thermo = thermodynamics;
    varset.reference_state = refstate;
    varset.primal_geometry = primal_geom;
    varset.dual_geometry = dual_geom;

    // If more physics parameterizations are added this logic might need to
    // change
    // varset.couple_wind = !(coupler.get_option<std::string>("sgs") == "none");
    varset.couple_wind = true;

    // dens_pos IS NOT BEING PROPERLY DEALLOCATED AT THE END OF THE RUN IE WHEN
    // THE POOL IS DESTROYED THIS IS REALLY WEIRD
    //  Allocate device arrays for whether densities are positive-definite
    // this->dens_pos       = bool1d("dens_pos"      ,ndensity);
    // boolHost1d dens_pos_host      ("dens_pos_host"      ,ndensity);

    for (int l = ndensity_dycore; l < ndensity_nophysics; l++) {
      // dens_pos_host(l) = params.dycore_tracerpos[l-ndensity_dycore];
      varset.dens_pos(l) = params.dycore_tracerpos[l - ndensity_dycore];
    }

    std::vector<std::string> tracer_names_loc = coupler.get_tracer_names();
    bool water_vapor_found = false;
    for (int tr = 0; tr < ntracers_physics; tr++) {
      bool found, positive, adds_mass;
      std::string desc;
      coupler.get_tracer_info(tracer_names_loc[tr], desc, found, positive,
                              adds_mass);
      varset.dens_name[tr + ndensity_nophysics] = tracer_names_loc[tr];
      varset.dens_desc[tr + ndensity_nophysics] = desc;
      // dens_pos_host      (tr+ndensity_nophysics) = positive ;
      varset.dens_pos(tr + ndensity_nophysics) = positive;
      // varset.dens_pos(tr + ndensity_nophysics) = false;
      if (tracer_names_loc[tr] == std::string("water_vapor")) {
        varset.dm_id_vap = tr;
        water_vapor_found = true;
      }
      if (tracer_names_loc[tr] == std::string("cloud_liquid") ||
          tracer_names_loc[tr] == std::string("cloud_water")) {
        varset.dm_id_liq = tr;
        varset.liquid_found = true;
      }
      if (tracer_names_loc[tr] == std::string("cloud_ice")) {
        varset.dm_id_ice = tr;
        varset.ice_found = true;
      }
    }
    if (ntracers_physics > 0) {
      if (!water_vapor_found) {
        endrun("ERROR: processed registered tracers, and water_vapor was not "
               "found");
      }
    }

    for (int i = ndensity_dycore; i < ndensity_nophysics; i++) {
      varset.dens_name[i] = "Tracer" + std::to_string(i - ndensity_dycore);
      varset.dens_desc[i] =
          "Dycore Tracer" + std::to_string(i - ndensity_dycore);
    }

    for (int i = 0; i < ndensity; i++) {
      serial_print("Density" + std::to_string(i) + " Name: " +
                       varset.dens_name[i] + " Desc: " + varset.dens_desc[i] +
                       " Pos: " + std::to_string(varset.dens_pos(i)),
                   params.masterproc);
    }

    // dens_pos_host      .deep_copy_to(dens_pos      );
    yakl::fence();
  }

  real YAKL_INLINE get_total_density(const real5d &densvar, int k, int j, int i,
                                     int ks, int js, int is, int n) const {};
  real YAKL_INLINE get_total_density(const real3d &densvar, int k, int ks,
                                     int n) const {};
  real YAKL_INLINE get_dry_density(const real5d &densvar, int k, int j, int i,
                                   int ks, int js, int is, int n) const {};
  real YAKL_INLINE get_entropic_var(const real5d &densvar, int k, int j, int i,
                                    int ks, int js, int is, int n) const {};
  real YAKL_INLINE get_entropic_var(const real3d &densvar, int k, int ks,
                                    int n) const {};
  void YAKL_INLINE set_density(real dens, real dry_dens, const real5d &densvar,
                               int k, int j, int i, int ks, int js, int is,
                               int n) const {};
  void YAKL_INLINE set_entropic_density(real entropic_var_density,
                                        const real5d &densvar, int k, int j,
                                        int i, int ks, int js, int is,
                                        int n) const {};
  real YAKL_INLINE get_alpha(const real5d &densvar, int k, int j, int i, int ks,
                             int js, int is, int n) const {};
  real YAKL_INLINE get_alpha(const real3d &densvar, int k, int ks,
                             int n) const {};
  real YAKL_INLINE get_qd(const real5d &densvar, int k, int j, int i, int ks,
                          int js, int is, int n) const {};
  real YAKL_INLINE get_qd(const real3d &densvar, int k, int ks, int n) const {};
  real YAKL_INLINE get_qv(const real5d &densvar, int k, int j, int i, int ks,
                          int js, int is, int n) const {};
  real YAKL_INLINE get_qv(const real3d &densvar, int k, int ks, int n) const {};
  real YAKL_INLINE get_ql(const real5d &densvar, int k, int j, int i, int ks,
                          int js, int is, int n) const {};
  real YAKL_INLINE get_qi(const real5d &densvar, int k, int j, int i, int ks,
                          int js, int is, int n) const {};
  real YAKL_INLINE _water_dens(const real5d &densvar, int k, int j, int i,
                               int ks, int js, int is, int n) const {};
  real YAKL_INLINE _water_dens(const real3d &densvar, int k, int ks,
                               int n) const {};

  void convert_dynamics_to_coupler_state(PamCoupler &coupler,
                                         const FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars);
  void convert_coupler_to_dynamics_state(PamCoupler &coupler,
                                         FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars);
};

template <class T>
void VariableSetBase<T>::convert_dynamics_to_coupler_state(
    PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars) {

  if constexpr (T::couple) {
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readwrite();

    real4d dm_dens_dry = dm.get<real, 4>("density_dry");
    real4d dm_uvel = dm.get<real, 4>("uvel");
    real4d dm_vvel = dm.get<real, 4>("vvel");
    real4d dm_wvel = dm.get<real, 4>("wvel");
    real4d dm_temp = dm.get<real, 4>("temp");

    pam::MultipleFields<ntracers_physics, real4d> dm_tracers;
    for (int tr = 0; tr < ntracers_physics; tr++) {
      auto trac = dm.get<real, 4>(dens_name[tr + ndensity_nophysics]);
      dm_tracers.add_field(trac);
    }

    if (couple_wind) {
      parallel_for(
          "Dynamics to Coupler State winds",
          SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                          dual_topology.n_cells_x, dual_topology.nens),
          YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
            // IN 3D THIS IS MORE COMPLICATED
            real uvel_l =
                prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) /
                primal_geometry.get_area_10entity(k + pks, j + pjs, i + pis, n);
            real uvel_r = prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs,
                                                          i + pis + 1, n) /
                          primal_geometry.get_area_10entity(k + pks, j + pjs,
                                                            i + pis + 1, n);
            real wvel_d = 0.0_fp;
            real wvel_u = 0.0_fp;
            if (k == 0) {
              wvel_u = 2.0_fp *
                       prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs,
                                                       i + pis, n) /
                       primal_geometry.get_area_01entity(k + pks, j + pjs,
                                                         i + pis, n);
              wvel_d = 0.0_fp;
            } else if (k == (dual_topology.nl)) {
              wvel_u = 0.0_fp;
              wvel_d = 2.0_fp *
                       prog_vars.fields_arr[WVAR].data(0, k + pks - 1, j + pjs,
                                                       i + pis, n) /
                       primal_geometry.get_area_01entity(k + pks - 1, j + pjs,
                                                         i + pis, n);
            } else {
              wvel_u = prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs,
                                                       i + pis, n) /
                       primal_geometry.get_area_01entity(k + pks, j + pjs,
                                                         i + pis, n);
              wvel_d = prog_vars.fields_arr[WVAR].data(0, k + pks - 1, j + pjs,
                                                       i + pis, n) /
                       primal_geometry.get_area_01entity(k + pks - 1, j + pjs,
                                                         i + pis, n);
            }
            // EVENTUALLY FIX THIS FOR 3D...
            real vvel = 0.0_fp;

            dm_uvel(k, j, i, n) = (uvel_l + uvel_r) * 0.5_fp;
            dm_vvel(k, j, i, n) = vvel;
            dm_wvel(k, j, i, n) = (wvel_u + wvel_d) * 0.5_fp;
          });
    }

    // std::cout << "Coupler state" << std::endl;
    parallel_for(
        "Dynamics to Coupler State densities",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          real qd = get_qd(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks,
                           djs, dis, n);
          real qv = get_qv(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks,
                           djs, dis, n);
          real alpha = get_alpha(prog_vars.fields_arr[DENSVAR].data, k, j, i,
                                 dks, djs, dis, n);
          real entropic_var = get_entropic_var(
              prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);

          real ql = 0.0_fp;
          if (liquid_found) {
            ql = get_ql(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs,
                        dis, n);
          }

          real qi = 0.0_fp;
          if (ice_found) {
            qi = get_qi(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs,
                        dis, n);
          }

          real temp = thermo.compute_T(alpha, entropic_var, qd, qv, ql, qi);

          dm_dens_dry(k, j, i, n) =
              get_dry_density(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks,
                              djs, dis, n) /
              dual_geometry.get_area_11entity(k + dks, j + djs, i + dis, n);
          dm_temp(k, j, i, n) = temp;
          for (int tr = ndensity_nophysics; tr < ndensity; tr++) {
            dm_tracers(tr - ndensity_nophysics, k, j, i, n) =
                prog_vars.fields_arr[DENSVAR].data(tr, k + dks, j + djs,
                                                   i + dis, n) /
                dual_geometry.get_area_11entity(k + dks, j + djs, i + dis, n);
          }
        });
  }
}

template <class T>
void VariableSetBase<T>::convert_coupler_to_dynamics_state(
    PamCoupler &coupler, FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars) {

  if constexpr (T::couple) {
    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readonly();

    auto dm_dens_dry = dm.get<real const, 4>("density_dry");
    auto dm_uvel = dm.get<real const, 4>("uvel");
    auto dm_vvel = dm.get<real const, 4>("vvel");
    auto dm_wvel = dm.get<real const, 4>("wvel");
    auto dm_temp = dm.get<real const, 4>("temp");

    pam::MultipleFields<ntracers_physics, realConst4d> dm_tracers;
    for (int tr = 0; tr < ntracers_physics; tr++) {
      auto trac = dm.get<real const, 4>(dens_name[tr + ndensity_nophysics]);
      dm_tracers.add_field(trac);
    }

    parallel_for(
        "Coupler to Dynamics State Densities",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
          real temp = dm_temp(k, j, i, n);

          real dens_dry = dm_dens_dry(k, j, i, n);
          real dens_vap = dm_tracers(dm_id_vap, k, j, i, n);
          real dens_liq = 0.0_fp;
          real dens_ice = 0.0_fp;
          if (liquid_found) {
            dens_liq = dm_tracers(dm_id_liq, k, j, i, n);
          }
          if (ice_found) {
            dens_ice = dm_tracers(dm_id_ice, k, j, i, n);
          }
          real dens = dens_dry + dens_ice + dens_liq + dens_vap;

          real qd = dens_dry / dens;
          real qv = dens_vap / dens;
          real ql = dens_liq / dens;
          real qi = dens_ice / dens;

          real alpha = 1.0_fp / dens;
          real entropic_var =
              thermo.compute_entropic_var_from_T(alpha, temp, qd, qv, ql, qi);

#if !defined _AN && !defined _MAN
          set_density(dens * dual_geometry.get_area_11entity(k + dks, j + djs,
                                                             i + dis, n),
                      dens_dry * dual_geometry.get_area_11entity(
                                     k + dks, j + djs, i + dis, n),
                      prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs,
                      dis, n);
#endif
          set_entropic_density(
              entropic_var * dens *
                  dual_geometry.get_area_11entity(k + dks, j + djs, i + dis, n),
              prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);

          for (int tr = ndensity_nophysics; tr < ndensity; tr++) {
            prog_vars.fields_arr[DENSVAR].data(tr, k + dks, j + djs, i + dis,
                                               n) =
                dm_tracers(tr - ndensity_nophysics, k, j, i, n) *
                dual_geometry.get_area_11entity(k + dks, j + djs, i + dis, n);
          }
        });

    if (couple_wind) {
      if (couple_wind_exact_inverse) {
        parallel_for(
            "Coupler to Dynamics State Primal U",
            SimpleBounds<3>(primal_topology.ni, primal_topology.n_cells_y,
                            primal_topology.nens),
            YAKL_CLASS_LAMBDA(int k, int j, int n) {
              real x0 = 0;
              for (int i = 0; i < primal_topology.n_cells_x; ++i) {
                x0 += (i % 2 == 0 ? 1 : -1) * dm_uvel(k, j, i, n);
              }
              prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, pis, n) = x0;
              prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, pis, n) *=
                  primal_geometry.get_area_10entity(k + pks, j + pjs, pis, n);

              for (int i = 1; i < primal_topology.n_cells_x; ++i) {
                x0 = 2 * dm_uvel(k, j, i - 1, n) - x0;
                prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) = x0;
                prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) *=
                    primal_geometry.get_area_10entity(k + pks, j + pjs, i + pis,
                                                      n);
              }
            });
        parallel_for(
            "Coupler to Dynamics State Primal W",
            SimpleBounds<3>(primal_topology.n_cells_y,
                            primal_topology.n_cells_x, primal_topology.nens),
            YAKL_CLASS_LAMBDA(int j, int i, int n) {
              real x0 = dm_wvel(0, j, i, n);
              prog_vars.fields_arr[WVAR].data(0, pks, j + pjs, i + pis, n) = x0;
              prog_vars.fields_arr[WVAR].data(0, pks, j + pjs, i + pis, n) *=
                  primal_geometry.get_area_01entity(pks, j + pjs, i + pis, n);

              for (int k = 1; k < primal_topology.nl; ++k) {
                x0 = 2 * dm_wvel(k, j, i, n) - x0;
                prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) = x0;
                prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) *=
                    primal_geometry.get_area_01entity(k + pks, j + pjs, i + pis,
                                                      n);
              }
            });
      } else {
        parallel_for(
            "Coupler to Dynamics State Primal U",
            SimpleBounds<4>(primal_topology.ni, primal_topology.n_cells_y,
                            primal_topology.n_cells_x, primal_topology.nens),
            YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
              // periodic wrapping
              int il = i - 1;
              if (i == 0) {
                il = primal_topology.n_cells_x - 1;
              }
              prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, i + pis, n) =
                  (dm_uvel(k, j, il, n) + dm_uvel(k, j, i, n)) * 0.5_fp *
                  primal_geometry.get_area_10entity(k + pks, j + pjs, i + pis,
                                                    n);
            });

        // EVENTUALLY THIS NEEDS TO HAVE A FLAG ON IT!
        parallel_for(
            "Coupler to Dynamics State Primal W",
            SimpleBounds<4>(primal_topology.nl, primal_topology.n_cells_y,
                            primal_topology.n_cells_x, primal_topology.nens),
            YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
              prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs, i + pis, n) =
                  (dm_wvel(k, j, i, n) + dm_wvel(k + 1, j, i, n)) * 0.5_fp *
                  primal_geometry.get_area_01entity(k + pks, j + pjs, i + pis,
                                                    n);
            });
      }
    }
  }
}

#ifdef _SWE
template <>
void VariableSetBase<VS_SWE>::initialize(PamCoupler &coupler,
                                         ModelParameters &params,
                                         const ThermoPotential &thermodynamics,
                                         const ReferenceState &refstate,
                                         const Geometry<Straight> &primal_geom,
                                         const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "h";
  dens_desc[0] = "fluid height";
  dens_pos(0) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}
#endif

#ifdef _TSWE
template <>
void VariableSetBase<VS_TSWE>::initialize(PamCoupler &coupler,
                                          ModelParameters &params,
                                          const ThermoPotential &thermodynamics,
                                          const ReferenceState &refstate,
                                          const Geometry<Straight> &primal_geom,
                                          const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "h";
  dens_name[1] = "S";
  dens_desc[0] = "fluid height";
  dens_desc[1] = "bouyancy density";
  dens_pos(0) = false;
  dens_pos(1) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}
#endif

#ifdef _CE
template <>
void VariableSetBase<VS_CE>::initialize(PamCoupler &coupler,
                                        ModelParameters &params,
                                        const ThermoPotential &thermodynamics,
                                        const ReferenceState &refstate,
                                        const Geometry<Straight> &primal_geom,
                                        const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "rho";
  dens_name[1] = "S";
  dens_desc[0] = "fluid density";
  dens_desc[1] = "entropic variable density";
  dens_pos(0) = false;
  dens_pos(1) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}

template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_entropic_var(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  return densvar(1, k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_entropic_var(const real3d &densvar,
                                                          int k, int ks,
                                                          int n) const {
  return densvar(1, k + ks, n) / densvar(0, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_alpha(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {
  return dual_geometry.get_area_11entity(k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_alpha(const real3d &densvar, int k,
                                                   int ks, int n) const {
  return dual_geometry.get_area_11entity(k + ks, 0, 0, n) /
         densvar(0, k + ks, n);
}
#endif

#ifdef _AN
template <>
void VariableSetBase<VS_AN>::initialize(PamCoupler &coupler,
                                        ModelParameters &params,
                                        const ThermoPotential &thermodynamics,
                                        const ReferenceState &refstate,
                                        const Geometry<Straight> &primal_geom,
                                        const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "S";
  dens_desc[0] = "entropic variable density";
  dens_pos(0) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return reference_state.dens.data(MASSDENSINDX, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_total_density(
    const real3d &densvar, int k, int ks, int n) const {
  return densvar(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_entropic_var(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  return densvar(ENTROPICDENSINDX, k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_entropic_var(const real3d &densvar,
                                                          int k, int ks,
                                                          int n) const {
  return densvar(ENTROPICDENSINDX, k + ks, n) /
         densvar(MASSDENSINDX, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_alpha(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {
  return dual_geometry.get_area_11entity(k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_alpha(const real3d &densvar, int k,
                                                   int ks, int n) const {
  return dual_geometry.get_area_11entity(k + ks, 0, 0, n) /
         densvar(MASSDENSINDX, k + ks, n);
}
#endif

// We rely on physics packages ie micro to provide water species- must at least
// have vapor and cloud liquid
#ifdef _MAN
template <>
void VariableSetBase<VS_MAN>::initialize(PamCoupler &coupler,
                                         ModelParameters &params,
                                         const ThermoPotential &thermodynamics,
                                         const ReferenceState &refstate,
                                         const Geometry<Straight> &primal_geom,
                                         const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "S";
  dens_desc[0] = "entropic variable density";
  dens_pos(0) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_total_density(
    const real3d &densvar, int k, int ks, int n) const {
  return reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_entropic_var(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(ENTROPICDENSINDX, k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_entropic_var(
    const real3d &densvar, int k, int ks, int n) const {
  return reference_state.dens.data(ENTROPICDENSINDX, k + ks, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_alpha(const real5d &densvar,
                                                    int k, int j, int i, int ks,
                                                    int js, int is,
                                                    int n) const {
  return dual_geometry.get_area_11entity(k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_alpha(const real3d &densvar,
                                                    int k, int ks,
                                                    int n) const {
  return dual_geometry.get_area_11entity(k + ks, 0, 0, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qv(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return densvar(dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qv(const real3d &densvar, int k,
                                                 int ks, int n) const {
  return densvar(dm_id_vap + ndensity_nophysics, k + ks, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_ql(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return densvar(dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qi(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return densvar(dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::_water_dens(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  real vap_dens =
      densvar(dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liquid_found) {
    liq_dens =
        densvar(dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n);
  }
  if (ice_found) {
    ice_dens =
        densvar(dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n);
  }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::_water_dens(const real3d &densvar,
                                                      int k, int ks,
                                                      int n) const {
  real vap_dens = densvar(dm_id_vap + ndensity_nophysics, k + ks, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liquid_found) {
    liq_dens = densvar(dm_id_liq + ndensity_nophysics, k + ks, n);
  }
  if (ice_found) {
    ice_dens = densvar(dm_id_ice + ndensity_nophysics, k + ks, n);
  }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_dry_density(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  return (reference_state.dens.data(MASSDENSINDX, k + ks, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n));
}
template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qd(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return (reference_state.dens.data(MASSDENSINDX, k + ks, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n)) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qd(const real3d &densvar, int k,
                                                 int ks, int n) const {
  return (reference_state.dens.data(MASSDENSINDX, k + ks, n) -
          _water_dens(densvar, k, ks, n)) /
         reference_state.dens.data(MASSDENSINDX, k + ks, n);
}

template <>
void YAKL_INLINE VariableSetBase<VS_MAN>::set_density(real dens, real dryden,
                                                      const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  reference_state.dens.data(MASSDENSINDX, k + ks, n) = dens;
}
template <>
void YAKL_INLINE VariableSetBase<VS_MAN>::set_entropic_density(
    real entropic_var_density, const real5d &densvar, int k, int j, int i,
    int ks, int js, int is, int n) const {
  densvar(ENTROPICDENSINDX, k + ks, j + js, i + is, n) = entropic_var_density;
}
#endif

#ifdef _MCErho
template <>
void VariableSetBase<VS_MCE_rho>::initialize(
    PamCoupler &coupler, ModelParameters &params,
    const ThermoPotential &thermodynamics, const ReferenceState &refstate,
    const Geometry<Straight> &primal_geom, const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "rho";
  dens_name[1] = "S";
  dens_desc[0] = "fluid density";
  dens_desc[1] = "entropic variable density";
  dens_pos(0) = false;
  dens_pos(1) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_total_density(
    const real3d &densvar, int k, int ks, int n) const {
  return densvar(0, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_entropic_var(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(1, k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_entropic_var(
    const real3d &densvar, int k, int ks, int n) const {
  return densvar(1, k + ks, n) / densvar(0, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_alpha(const real5d &densvar,
                                                        int k, int j, int i,
                                                        int ks, int js, int is,
                                                        int n) const {
  return dual_geometry.get_area_11entity(k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_alpha(const real3d &densvar,
                                                        int k, int ks,
                                                        int n) const {
  return dual_geometry.get_area_11entity(k + ks, 0, 0, n) /
         densvar(0, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qv(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return densvar(dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qv(const real3d &densvar,
                                                     int k, int ks,
                                                     int n) const {
  return densvar(dm_id_vap + ndensity_nophysics, k + ks, n) /
         densvar(0, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_ql(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return densvar(dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qi(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return densvar(dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n) /
         densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::_water_dens(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  real vap_dens =
      densvar(dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liquid_found) {
    liq_dens =
        densvar(dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n);
  }
  if (ice_found) {
    ice_dens =
        densvar(dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n);
  }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::_water_dens(const real3d &densvar,
                                                          int k, int ks,
                                                          int n) const {
  real vap_dens = densvar(dm_id_vap + ndensity_nophysics, k + ks, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liquid_found) {
    liq_dens = densvar(dm_id_liq + ndensity_nophysics, k + ks, n);
  }
  if (ice_found) {
    ice_dens = densvar(dm_id_ice + ndensity_nophysics, k + ks, n);
  }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_dry_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return (densvar(0, k + ks, j + js, i + is, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n));
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qd(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return (densvar(0, k + ks, j + js, i + is, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n)) /
         densvar(0, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qd(const real3d &densvar,
                                                     int k, int ks,
                                                     int n) const {
  return (densvar(0, k + ks, n) - _water_dens(densvar, k, ks, n)) /
         densvar(0, k + ks, n);
}

template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rho>::set_density(
    real dens, real dryden, const real5d &densvar, int k, int j, int i, int ks,
    int js, int is, int n) const {
  densvar(0, k + ks, j + js, i + is, n) = dens;
}
template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rho>::set_entropic_density(
    real entropic_var_density, const real5d &densvar, int k, int j, int i,
    int ks, int js, int is, int n) const {
  densvar(1, k + ks, j + js, i + is, n) = entropic_var_density;
}
#endif

#ifdef _MCErhod
template <>
void VariableSetBase<VS_MCE_rhod>::initialize(
    PamCoupler &coupler, ModelParameters &params,
    const ThermoPotential &thermodynamics, const ReferenceState &refstate,
    const Geometry<Straight> &primal_geom, const Geometry<Twisted> &dual_geom) {
  dens_name[0] = "rho_d";
  dens_name[1] = "S";
  dens_desc[0] = "fluid dry density";
  dens_desc[1] = "entropic variable density";
  dens_pos(0) = false;
  dens_pos(1) = false;
  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  real dry_dens = densvar(0, k + ks, j + js, i + is, n);
  real vap_dens =
      densvar(dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liquid_found) {
    liq_dens =
        densvar(dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n);
  }
  if (ice_found) {
    ice_dens =
        densvar(dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n);
  }
  return dry_dens + vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_dry_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(0, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_entropic_var(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(1, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_alpha(const real5d &densvar,
                                                         int k, int j, int i,
                                                         int ks, int js, int is,
                                                         int n) const {
  return dual_geometry.get_area_11entity(k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_qd(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(0, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_qv(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dm_id_vap + ndensity_nophysics, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_ql(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dm_id_liq + ndensity_nophysics, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_qi(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dm_id_ice + ndensity_nophysics, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rhod>::set_density(
    real dens, real drydens, const real5d &densvar, int k, int j, int i, int ks,
    int js, int is, int n) const {
  densvar(0, k + ks, j + js, i + is, n) = drydens;
}
template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rhod>::set_entropic_density(
    real entropic_var_density, const real5d &densvar, int k, int j, int i,
    int ks, int js, int is, int n) const {
  densvar(1, k + ks, j + js, i + is, n) = entropic_var_density;
}
#endif

#ifdef _SWE
using VariableSet = VariableSetBase<VS_SWE>;
#elif _TSWE
using VariableSet = VariableSetBase<VS_TSWE>;
#elif _CE
using VariableSet = VariableSetBase<VS_CE>;
#elif _AN
using VariableSet = VariableSetBase<VS_AN>;
#elif _MAN
using VariableSet = VariableSetBase<VS_MAN>;
#elif _MCErho
using VariableSet = VariableSetBase<VS_MCE_rho>;
#elif _MCErhod
using VariableSet = VariableSetBase<VS_MCE_rhod>;
#elif _CEp
using VariableSet = VariableSetBase<VS_CE_p>;
#elif _MCErhop
using VariableSet = VariableSetBase<VS_MCE_rhop>;
#elif _MCErhodp
using VariableSet = VariableSetBase<VS_MCE_rhodp>;
#endif
