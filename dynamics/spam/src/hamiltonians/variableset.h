#pragma once

#include "DataManager.h"
#include "Microphysics.h"
#include "MultipleFields.h"
#include "SGS.h"
#include "common.h"
#include "geometry.h"
#include "pam_coupler.h" //Has DataManager and pam_const
#include "refstate.h"
#include "thermo.h"
using pam::PamCoupler;

//#ifdef PAM_STANDALONE
  // solve a system to exactly invert the velocity averaging done
  // during conversion to coupler state when coupling winds
//  constexpr bool couple_wind_exact_inverse = false;
//#else
  // disable exact inversion of velocity averaging
  // for more flexible MMF configuration
//  constexpr bool couple_wind_exact_inverse = false;
//#endif

struct VS_SWE {
  static constexpr bool couple = false;

  static constexpr uint ndensity_dycore = 1;
  static constexpr uint ndensity_dycore_prognostic = ndensity_dycore;
  static constexpr uint ndensity_dycore_active = ndensity_dycore;

  static constexpr uint ntracers_dycore_active =
      std::min<uint>(3, ntracers_dycore);

  static constexpr uint ntracers_physics = 0;
  static constexpr uint ntracers_physics_active = 0;
};

struct VS_TSWE {
  static constexpr bool couple = false;

  static constexpr uint ndensity_dycore = 2;
  static constexpr uint ndensity_dycore_prognostic = 2;
  static constexpr uint ndensity_dycore_active = 2;

  static constexpr uint ntracers_dycore_active = 0;

  static constexpr uint ntracers_physics = 0;
  static constexpr uint ntracers_physics_active = 0;
};

struct VS_CE {
  static constexpr bool couple = false;

  static constexpr uint ndensity_dycore = 2;
  static constexpr uint ndensity_dycore_prognostic = ndensity_dycore;
  static constexpr uint ndensity_dycore_active = ndensity_dycore;

  static constexpr uint ntracers_dycore_active = 0;

  static constexpr uint ntracers_physics = 0;
  static constexpr uint ntracers_physics_active = 0;
};

struct VS_AN {
  static constexpr bool couple = false;

  static constexpr uint ndensity_dycore = 2;
  static constexpr uint ndensity_dycore_prognostic = 1;
  static constexpr uint ndensity_dycore_active = ndensity_dycore;

  static constexpr uint ntracers_dycore_active = 0;

  static constexpr uint ntracers_physics = 0;
  static constexpr uint ntracers_physics_active = 0;
};

struct VS_MAN {
  static constexpr bool couple = true;

  static constexpr uint nmoist = 4;

  static constexpr uint ndensity_dycore = 2;
  static constexpr uint ndensity_dycore_prognostic = 1;
  static constexpr uint ndensity_dycore_active = ndensity_dycore;

  static constexpr uint ntracers_dycore_active = 0;

  static constexpr uint ntracers_physics =
      Microphysics::get_num_tracers() + SGS::get_num_tracers();
  static constexpr uint ntracers_physics_active =
      ThermoPotential::moist_species_decouple_from_dynamics
          ? 0
          : std::min<uint>(3, Microphysics::get_num_tracers());
};

struct VS_MCE_rho {
  static constexpr bool couple = true;

  static constexpr uint nmoist = 4;

  static constexpr uint ndensity_dycore = 2;
  static constexpr uint ndensity_dycore_prognostic = ndensity_dycore;
  static constexpr uint ndensity_dycore_active = ndensity_dycore;

  static constexpr uint ntracers_dycore_active = 0;

  static constexpr uint ntracers_physics =
      Microphysics::get_num_tracers() + SGS::get_num_tracers();
  static constexpr uint ntracers_physics_active =
      ThermoPotential::moist_species_decouple_from_dynamics
          ? 0
          : std::min<uint>(3, Microphysics::get_num_tracers());
};

struct VS_MCE_rhod : VS_MCE_rho {};

struct VS_CE_p {
  static constexpr bool couple = false;
};

struct VS_MCE_rhop {
  static constexpr bool couple = false;
};

struct VS_MCE_rhodp {
  static constexpr bool couple = false;
};

template <class T> class VariableSetBase : public T {
public:
  using T::ndensity_dycore;
  using T::ndensity_dycore_active;
  using T::ndensity_dycore_prognostic;

  using T::ntracers_dycore_active;
  static_assert(ntracers_dycore_active <= ntracers_dycore);

  using T::ntracers_physics;
  using T::ntracers_physics_active;

  static constexpr uint ndensity =
      ndensity_dycore + ntracers_dycore + ntracers_physics;
  static constexpr uint ndensity_nophysics =
      ndensity_dycore_prognostic + ntracers_dycore;
  static constexpr uint ndensity_active =
      ndensity_dycore_active + ntracers_dycore_active + ntracers_physics_active;
  static constexpr uint ndensity_prognostic =
      ndensity_dycore_prognostic + ntracers_dycore + ntracers_physics;

  std::string dens_name[ndensity]; // Name of each density
  std::string dens_desc[ndensity]; // Description of each density
  SArray<bool, 1, ndensity_prognostic>
      dens_pos; // Whether each density is positive-definite
  SArray<bool, 1, ndensity> dens_active; // Whether each density is active
  SArray<bool, 1, ndensity>
      dens_prognostic; // Whether each density is prognostic
  SArray<int, 1, ndensity_active>
      active_dens_ids; // indices of active densities
  // bool couple_wind;

  int dm_id_vap = std::numeric_limits<int>::min();
  int dm_id_liq = std::numeric_limits<int>::min();
  int dm_id_ice = std::numeric_limits<int>::min();

  int dens_id_mass = std::numeric_limits<int>::min();
  int dens_id_entr = std::numeric_limits<int>::min();
  int dens_id_vap = std::numeric_limits<int>::min();
  int dens_id_liq = std::numeric_limits<int>::min();
  int dens_id_ice = std::numeric_limits<int>::min();

  int active_id_mass = std::numeric_limits<int>::min();
  int active_id_entr = std::numeric_limits<int>::min();
  int active_id_vap = std::numeric_limits<int>::min();
  int active_id_liq = std::numeric_limits<int>::min();
  int active_id_ice = std::numeric_limits<int>::min();

  bool ice_found = false;
  bool liq_found = false;

  ThermoPotential thermo;
  ReferenceState reference_state;
  Geometry<Straight> primal_geometry;
  Geometry<Twisted> dual_geometry;

  void initialize(PamCoupler &coupler, ModelParameters &params,
                  const ThermoPotential &thermodynamics,
                  const ReferenceState &refstate,
                  const Geometry<Straight> &primal_geom,
                  const Geometry<Twisted> &dual_geom,
                  bool verbose=false);
  static void initialize(VariableSetBase &varset, PamCoupler &coupler,
                         ModelParameters &params,
                         const ThermoPotential &thermodynamics,
                         const ReferenceState &refstate,
                         const Geometry<Straight> &primal_geom,
                         const Geometry<Twisted> &dual_geom,
                         bool verbose=false) {

    if (T::couple && params.couple_wind_exact_inverse) {
      if (primal_geom.topology.n_cells_x % 2 == 0) {
        throw std::runtime_error(
            "The number of crm cells in the horizontal "
            "has to be odd when using the couple_wind_exact_inverse option");
      }
    }
    varset.thermo = thermodynamics;
    varset.reference_state = refstate;
    varset.primal_geometry = primal_geom;
    varset.dual_geometry = dual_geom;

    // Need different logic for this.
    // maybe a coupler option so that the driver can set it directly??
    //varset.couple_wind = !(coupler.get_option<std::string>("sgs") == "none") ||
    //                     !(coupler.option_exists("standalone_input_file"));
    varset.couple_wind = true;

    for (int l = ndensity_dycore_prognostic; l < ndensity_nophysics; l++) {
      varset.dens_pos(l) =
          params.dycore_tracerpos[l - ndensity_dycore_prognostic];
      varset.dens_prognostic(l) = true;
      varset.dens_active(l) =
          (l - ndensity_dycore_prognostic) < ntracers_dycore_active ? true
                                                                    : false;
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
      varset.dens_pos(tr + ndensity_nophysics) = positive;
      varset.dens_prognostic(tr + ndensity_nophysics) = true;
      varset.dens_active(tr + ndensity_nophysics) = false;
      if (tracer_names_loc[tr] == std::string("water_vapor")) {
        varset.dm_id_vap = tr;
        varset.dens_id_vap = ndensity_nophysics + tr;
        water_vapor_found = true;
        if (!::ThermoPotential::moist_species_decouple_from_dynamics) {
          varset.dens_active(tr + ndensity_nophysics) = true;
        }
      }
      if (tracer_names_loc[tr] == std::string("cloud_liquid") ||
          tracer_names_loc[tr] == std::string("cloud_water")) {
        varset.dm_id_liq = tr;
        varset.dens_id_liq = ndensity_nophysics + tr;
        varset.liq_found = true;
        if (!::ThermoPotential::moist_species_decouple_from_dynamics) {
          varset.dens_active(tr + ndensity_nophysics) = true;
        }
      }
      if (tracer_names_loc[tr] == std::string("ice")) {
        varset.dm_id_ice = tr;
        varset.dens_id_ice = ndensity_nophysics + tr;
        varset.ice_found = true;
        if (!::ThermoPotential::moist_species_decouple_from_dynamics) {
          varset.dens_active(tr + ndensity_nophysics) = true;
        }
      }
    }
    if (ntracers_physics > 0) {
      if (!water_vapor_found) {
        endrun("ERROR: processed registered tracers, and water_vapor was not "
               "found");
      }
    }

    // get indicies of active densities
    int active_i = 0;
    for (int i = 0; i < ndensity; ++i) {
      if (varset.dens_active(i)) {
        varset.active_dens_ids(active_i) = i;
        if (water_vapor_found && i == ndensity_nophysics + varset.dm_id_vap) {
          varset.active_id_vap = active_i;
        }
        if (varset.ice_found && i == ndensity_nophysics + varset.dm_id_ice) {
          varset.active_id_ice = active_i;
        }
        if (varset.liq_found && i == ndensity_nophysics + varset.dm_id_liq) {
          varset.active_id_liq = active_i;
        }
        active_i++;
      }
    }

    for (int i = ndensity_dycore_prognostic; i < ndensity_nophysics; i++) {
      varset.dens_name[i] =
          "Tracer" + std::to_string(i - ndensity_dycore_prognostic);
      varset.dens_desc[i] =
          "Dycore Tracer" + std::to_string(i - ndensity_dycore_prognostic);
    }

    if (verbose) {
      serial_print("PAM-C densities", params.masterproc);
      for (int i = 0; i < ndensity; i++) {
        std::stringstream ss;
        ss << std::left;
        ss << std::setw(4) << i;
        ss << std::setw(21) << varset.dens_name[i].substr(0, 19);
        ss << std::setw(29) << varset.dens_desc[i].substr(0, 27);
        ss << std::setw(12) << (varset.dens_prognostic(i) ? "prognostic" : "");
        ss << std::setw(8) << (varset.dens_active(i) ? "active" : "");
        ss << std::setw(10)
           << ((varset.dens_prognostic(i) && varset.dens_pos(i)) ? "positive"
                                                                 : "");

        serial_print(ss.str(), params.masterproc);
      }
    }

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

  real YAKL_INLINE get_pres(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {};
  real YAKL_INLINE get_pres(int k, int ks, int n) const {};
  real YAKL_INLINE get_ref_dens(int k, int ks, int n) const {};

  void pamc_debug_chk(int id, PamCoupler &coupler,
                                         const FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nprognostic> &prev_vars);

  void convert_dynamics_to_coupler_densities(PamCoupler &coupler,
                                         const FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars);
  void convert_dynamics_to_coupler_wind(PamCoupler &coupler,
                                         const FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars, bool couple_wind_exact_inverse);
  void convert_dynamics_to_coupler_staggered_wind(PamCoupler &coupler,
                                         const FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars);
  void convert_coupler_to_dynamics_densities(PamCoupler &coupler,
                                         FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars);
  void convert_coupler_to_dynamics_wind(PamCoupler &coupler,
                                         FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars, bool couple_wind_exact_inverse);
  void convert_coupler_to_dynamics_staggered_wind(PamCoupler &coupler,
                                         FieldSet<nprognostic> &prog_vars,
                                         const FieldSet<nconstant> &const_vars);
};

//THIS IS ANELASTIC SPECIFIC FOR NOW, IDEALLY IT SHOULD MADE TO WORK FOR EITHER AN OR COMPRESSIBLE
template <class T>
void VariableSetBase<T>::pamc_debug_chk(int id,
                                        PamCoupler &coupler,
                                        const FieldSet<nprognostic> &prog_vars,
                                        const FieldSet<nprognostic> &prev_vars)
{
  const auto &primal_topology = primal_geometry.topology;
  const auto &dual_topology = dual_geometry.topology;
  int dis = dual_topology.is;
  int djs = dual_topology.js;
  int dks = dual_topology.ks;
  int pis = primal_topology.is;
  int pjs = primal_topology.js;
  int pks = primal_topology.ks;
  // const auto densvar = prog_vars.fields_arr[DENSVAR].data;
  parallel_for("pamc_debug_chk",SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y, dual_topology.n_cells_x, dual_topology.nens), YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {

    const auto densvar1 = prev_vars.fields_arr[DENSVAR].data;
    real entr1 = get_entropic_var( densvar1, k, j, i, dks, djs, dis, n);
    real pres1 = get_pres(         densvar1, k, j, i, dks, djs, dis, n);
    real qv1   = get_qv(           densvar1, k, j, i, dks, djs, dis, n);
    real ql1   = get_ql(           densvar1, k, j, i, dks, djs, dis, n);
    real qi1   = get_qi(           densvar1, k, j, i, dks, djs, dis, n);
    real qd1   = 1.0_fp - qv1 - ql1 - qi1;
    real temp1 = thermo.compute_T_from_p(pres1, entr1, qd1, qv1, ql1, qi1);

    const auto densvar2 = prog_vars.fields_arr[DENSVAR].data;
    real entr2 = get_entropic_var( densvar2, k, j, i, dks, djs, dis, n);
    real pres2 = get_pres(         densvar2, k, j, i, dks, djs, dis, n);
    real qv2   = get_qv(           densvar2, k, j, i, dks, djs, dis, n);
    real ql2   = get_ql(           densvar2, k, j, i, dks, djs, dis, n);
    real qi2   = get_qi(           densvar2, k, j, i, dks, djs, dis, n);
    real qd2   = 1.0_fp - qv2 - ql2 - qi2;
    real temp2 = thermo.compute_T_from_p(pres2, entr2, qd2, qv2, ql2, qi2);

    if ( isnan(temp2) || isnan(qv2) || isnan(ql2) || isnan(qi2) || isnan(qd2) || temp2<0 || entr2<0 ) {
      printf("pamc_debug_chk - id:%d k:%d j:%d i:%d n:%d en_var :: %g  =>  %g \n",id,k,j,i,n,entr1,entr2);
      printf("pamc_debug_chk - id:%d k:%d j:%d i:%d n:%d temp   :: %g  =>  %g \n",id,k,j,i,n,temp1,temp2);
      printf("pamc_debug_chk - id:%d k:%d j:%d i:%d n:%d qd     :: %g  =>  %g \n",id,k,j,i,n,qd1,qd2);
      printf("pamc_debug_chk - id:%d k:%d j:%d i:%d n:%d qv     :: %g  =>  %g \n",id,k,j,i,n,qv1,qv2);
      printf("pamc_debug_chk - id:%d k:%d j:%d i:%d n:%d ql     :: %g  =>  %g \n",id,k,j,i,n,ql1,ql2);
      printf("pamc_debug_chk - id:%d k:%d j:%d i:%d n:%d qi     :: %g  =>  %g \n",id,k,j,i,n,qi1,qi2);
    }
  });
}


template <class T>
void VariableSetBase<T>::convert_dynamics_to_coupler_densities(
    PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars) {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    auto &dm = coupler.get_data_manager_device_readwrite();

  if constexpr (T::couple) {

    real4d dm_dens_dry = dm.get<real, 4>("density_dry");
    real4d dm_temp = dm.get<real, 4>("temp");

    pam::MultipleFields<ntracers_physics, real4d> dm_tracers;
    for (int tr = 0; tr < ntracers_physics; tr++) {
      auto trac = dm.get<real, 4>(dens_name[tr + ndensity_nophysics]);
      dm_tracers.add_field(trac);
    }

    parallel_for(
        "Dynamics to Coupler State densities",
        SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                        dual_topology.n_cells_x, dual_topology.nens),
        YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {


          real qd = get_qd(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks,
                          djs, dis, n);
          real qv = get_qv(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks,
                           djs, dis, n);

          real entropic_var = get_entropic_var(
              prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);

          real ql = 0.0_fp;
          if (liq_found) {
            ql = get_ql(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs,
                        dis, n);
          }

          real qi = 0.0_fp;
          if (ice_found) {
            qi = get_qi(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs,
                        dis, n);
          }


//THIS LOGIC SHOULD PROBABLY LIVE IN VARIABLE SET, BUT IT IS NOT CLEAR HOW TO CLEANLY INTEGRATE IT YET
//Something like:
//set_dry_dens(dm_dry_dens,DENSARR,k,j,i,dks,djs,dis,n)
//compute_temp(DENSARR,k,j,i,dks,djs,dis,n, entropic_var, qd, qv, ql, qi) or compute_temp(total_dens, entropic_var, qd, qv, ql, qi) and get_total_dens(DENSARR,...)

#if defined(_CE) || defined(_MCErho) || defined(_MCErhod) || defined(_CEp) || defined(_MCErhop) || defined(_MCErhodp)
          dm_dens_dry(k, j, i, n) =
              get_dry_density(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks,
                              djs, dis, n) /
              dual_geometry.get_area_n1entity(k + dks, j + djs, i + dis, n);
#endif

#if defined(_CE) || defined(_MCErho) || defined(_MCErhod)
          real alpha = get_alpha(prog_vars.fields_arr[DENSVAR].data, k, j, i,
                                 dks, djs, dis, n);
          real temp = thermo.compute_T_from_alpha(alpha, entropic_var, qd, qv, ql, qi);
#elif defined(_CEp) || defined(_MCErhop) || defined(_MCErhodp) || defined(_AN) || defined(_MAN)
          real p = get_pres(prog_vars.fields_arr[DENSVAR].data, k, j, i,
                                 dks, djs, dis, n);
          real temp = thermo.compute_T_from_p(p, entropic_var, qd, qv, ql, qi);
#endif

          dm_temp(k, j, i, n) = temp;

          for (int tr = ndensity_nophysics; tr < ndensity_prognostic; tr++) {
            dm_tracers(tr - ndensity_nophysics, k, j, i, n) =
                prog_vars.fields_arr[DENSVAR].data(tr, k + dks, j + djs,
                                                   i + dis, n) /
                dual_geometry.get_area_n1entity(k + dks, j + djs, i + dis, n);
          }

        });

}
}

template <class T>
void VariableSetBase<T>::convert_dynamics_to_coupler_wind(
    PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars, bool couple_wind_exact_inverse) {

    const auto &primal_topology = primal_geometry.topology;
    const auto &dual_topology = dual_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readwrite();

  if constexpr (T::couple) {

    real4d dm_uvel = dm.get<real, 4>("uvel");
    real4d dm_vvel = dm.get<real, 4>("vvel");
    real4d dm_wvel = dm.get<real, 4>("wvel");

      parallel_for(
          "Dynamics to Coupler State winds",
          SimpleBounds<4>(dual_topology.nl, dual_topology.n_cells_y,
                          dual_topology.n_cells_x, dual_topology.nens),
          YAKL_CLASS_LAMBDA(int k, int j, int i, int n) {
            // IN 3D THIS IS MORE COMPLICATED
            real uvel_l = prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs,
                                                          i + pis, n) /
                          primal_geometry.get_area_10entity(0, k + pks, j + pjs,
                                                            i + pis, n);
            real uvel_r = prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs,
                                                          i + pis + 1, n) /
                          primal_geometry.get_area_10entity(0, k + pks, j + pjs,
                                                            i + pis + 1, n);
            real wvel_mid;
            if (k == 0) {
              wvel_mid = prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs,
                                                         i + pis, n) /
                         primal_geometry.get_area_01entity(k + pks, j + pjs,
                                                           i + pis, n);
            } else if (k == (dual_topology.nl)) {
              wvel_mid = prog_vars.fields_arr[WVAR].data(0, k + pks - 1,
                                                         j + pjs, i + pis, n) /
                         primal_geometry.get_area_01entity(k + pks - 1, j + pjs,
                                                           i + pis, n);
            } else {

              real e_u = primal_geometry.get_area_01entity(k + pks, j + pjs,
                                                           i + pis, n);
              real e_d = primal_geometry.get_area_01entity(k - 1 + pks, j + pjs,
                                                           i + pis, n);

              real wvel_u = prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs,
                                                            i + pis, n) /
                            e_u;
              real wvel_d = prog_vars.fields_arr[WVAR].data(
                                0, k + pks - 1, j + pjs, i + pis, n) /
                            e_d;

              wvel_mid = wvel_d + (wvel_u - wvel_d) * e_d / (e_u + e_d);
            }
            // EVENTUALLY FIX THIS FOR 3D...
            real vvel = 0.0_fp;

            dm_uvel(k, j, i, n) = (uvel_l + uvel_r) * 0.5_fp;
            dm_vvel(k, j, i, n) = vvel;
            dm_wvel(k, j, i, n) = wvel_mid;
          });
}
}


template <class T>
void VariableSetBase<T>::convert_dynamics_to_coupler_staggered_wind(
    PamCoupler &coupler, const FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars) {


    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readwrite();
  if constexpr (T::couple) {

//ADD THIS
}

}



template <class T>
void VariableSetBase<T>::convert_coupler_to_dynamics_densities(
    PamCoupler &coupler, FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars) {

  if constexpr (T::couple) {

    const auto &dual_topology = dual_geometry.topology;

    int dis = dual_topology.is;
    int djs = dual_topology.js;
    int dks = dual_topology.ks;

    auto &dm = coupler.get_data_manager_device_readonly();

    auto dm_dens_dry = dm.get<real const, 4>("density_dry");
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

          real dens_vap = dm_tracers(dm_id_vap, k, j, i, n);
          real dens_liq = 0.0_fp;
          real dens_ice = 0.0_fp;
          if (liq_found) {
            dens_liq = dm_tracers(dm_id_liq, k, j, i, n);
          }
          if (ice_found) {
            dens_ice = dm_tracers(dm_id_ice, k, j, i, n);
          }

//AGAIN, THIS LOGIC SHOULD REALLY LIVE IN VARIABLE SET
//Something like:
//get_total_dens(dens,dry_dens)
//get_qs(dens,dry_dens,dens_s)
//compute_entropic_var(dens,dry_dens,temp,qs)

#if defined(_CE) || defined(_MCErho) || defined(_MCE_rhod) || defined(_CEp)  || defined(_MCErhop) || defined(_MCErhodp)

          real dens_dry = dm_dens_dry(k, j, i, n);
          real dens = dens_dry + dens_vap;

          set_density(dens * dual_geometry.get_area_n1entity(k + dks, j + djs,
                                                             i + dis, n),
                      dens_dry * dual_geometry.get_area_n1entity(
                                     k + dks, j + djs, i + dis, n),
                      prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs,
                      dis, n);

          real qd = dens_dry / dens;
          real qv = dens_vap / dens;
          real ql = dens_liq / dens;
          real qi = dens_ice / dens;
          real alpha = 1.0_fp / dens;

#elif defined(_AN) || defined(_MAN)
          real dens = get_ref_dens(k, dks, n) / dual_geometry.get_area_n1entity(k + dks, j + djs, i + dis, n);

          real qv = dens_vap / dens;
          real ql = dens_liq / dens;
          real qi = dens_ice / dens;

          real qd = get_qd(prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);
#endif


#if defined(_CE) || defined(_MCErho) || defined(_MCE_rhod)
          real entropic_var =
              thermo.compute_entropic_var_from_T_alpha(alpha, temp, qd, qv, ql, qi);
#elif defined(_CEp)  || defined(_MCErhop) || defined(_MCErhodp)
          real p = get_pres(alpha, temp, qd, qv, ql, qi);
          real entropic_var =
              thermo.compute_entropic_var_from_T_p(p, temp, qd, qv, ql, qi);
#elif defined(_AN) || defined(_MAN)
          real p = get_pres(k, dks, n);
          real entropic_var =
              thermo.compute_entropic_var_from_T_p(p, temp, qd, qv, ql, qi);
#endif

          set_entropic_density(
              entropic_var * dens *
                  dual_geometry.get_area_n1entity(k + dks, j + djs, i + dis, n),
              prog_vars.fields_arr[DENSVAR].data, k, j, i, dks, djs, dis, n);

          for (int tr = ndensity_nophysics; tr < ndensity_prognostic; tr++) {
            prog_vars.fields_arr[DENSVAR].data(tr, k + dks, j + djs, i + dis,
                                               n) =
                dm_tracers(tr - ndensity_nophysics, k, j, i, n) *
                dual_geometry.get_area_n1entity(k + dks, j + djs, i + dis, n);
          }
        });
  }
}


template <class T>
void VariableSetBase<T>::convert_coupler_to_dynamics_wind(
    PamCoupler &coupler, FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars, bool couple_wind_exact_inverse) {

  if constexpr (T::couple) {

    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readonly();

    auto dm_uvel = dm.get<real const, 4>("uvel");
    auto dm_vvel = dm.get<real const, 4>("vvel");
    auto dm_wvel = dm.get<real const, 4>("wvel");

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
                  primal_geometry.get_area_10entity(0, k + pks, j + pjs, pis,
                                                    n);

              for (int i = 1; i < primal_topology.n_cells_x; ++i) {
                x0 = 2 * dm_uvel(k, j, i - 1, n) - x0;
                prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) = x0;
                prog_vars.fields_arr[VVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) *=
                    primal_geometry.get_area_10entity(0, k + pks, j + pjs,
                                                      i + pis, n);
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

                real ek = primal_geometry.get_area_01entity(k + pks, j + pjs,
                                                            i + pis, n);
                real ekm1 = primal_geometry.get_area_01entity(
                    k - 1 + pks, j + pjs, i + pis, n);
                x0 = (ek + ekm1) / ekm1 * dm_wvel(k, j, i, n) - x0 * ek / ekm1;
                prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) = x0;
                prog_vars.fields_arr[WVAR].data(0, k + pks, j + pjs, i + pis,
                                                n) *= ek;
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
                  primal_geometry.get_area_10entity(0, k + pks, j + pjs,
                                                    i + pis, n);
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


template <class T>
void VariableSetBase<T>::convert_coupler_to_dynamics_staggered_wind(
    PamCoupler &coupler, FieldSet<nprognostic> &prog_vars,
    const FieldSet<nconstant> &const_vars) {

  if constexpr (T::couple) {
    const auto &primal_topology = primal_geometry.topology;

    int pis = primal_topology.is;
    int pjs = primal_topology.js;
    int pks = primal_topology.ks;

    auto &dm = coupler.get_data_manager_device_readonly();

    auto dm_uvel = dm.get<real const, 4>("uvel");
    auto dm_vvel = dm.get<real const, 4>("vvel");
    auto dm_wvel = dm.get<real const, 4>("wvel");

//ADD THIS

      }


  }


#ifdef _SWE
template <>
void VariableSetBase<VS_SWE>::initialize(PamCoupler &coupler,
                                         ModelParameters &params,
                                         const ThermoPotential &thermodynamics,
                                         const ReferenceState &refstate,
                                         const Geometry<Straight> &primal_geom,
                                         const Geometry<Twisted> &dual_geom,
                                         bool verbose) {
  dens_id_mass = 0;
  active_id_mass = 0;
  dens_name[dens_id_mass] = "h";
  dens_desc[dens_id_mass] = "fluid height";
  dens_prognostic(dens_id_mass) = true;
  dens_active(dens_id_mass) = true;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}

template <>
real YAKL_INLINE VariableSetBase<VS_SWE>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_mass, k + ks, j + js, i + is, n);
}
#endif

#ifdef _TSWE
template <>
void VariableSetBase<VS_TSWE>::initialize(PamCoupler &coupler,
                                          ModelParameters &params,
                                          const ThermoPotential &thermodynamics,
                                          const ReferenceState &refstate,
                                          const Geometry<Straight> &primal_geom,
                                          const Geometry<Twisted> &dual_geom,
                                          bool verbose) {
  dens_id_mass = 0;
  active_id_mass = 0;
  dens_name[dens_id_mass] = "h";
  dens_desc[dens_id_mass] = "fluid height";
  dens_prognostic(dens_id_mass) = true;
  dens_active(dens_id_mass) = true;
  dens_pos(dens_id_mass) = false;

  dens_id_entr = 1;
  active_id_entr = 1;
  dens_name[dens_id_entr] = "S";
  dens_desc[dens_id_entr] = "bouyancy density";
  dens_prognostic(dens_id_entr) = true;
  dens_active(dens_id_entr) = true;
  dens_pos(dens_id_entr) = false;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}
#endif

#ifdef _CE
template <>
void VariableSetBase<VS_CE>::initialize(PamCoupler &coupler,
                                        ModelParameters &params,
                                        const ThermoPotential &thermodynamics,
                                        const ReferenceState &refstate,
                                        const Geometry<Straight> &primal_geom,
                                        const Geometry<Twisted> &dual_geom,
                                        bool verbose) {

  dens_id_mass = 0;
  active_id_mass = 0;
  dens_name[dens_id_mass] = "rho";
  dens_desc[dens_id_mass] = "fluid density";
  dens_prognostic(dens_id_mass) = true;
  dens_active(dens_id_mass) = true;
  dens_pos(dens_id_mass) = false;

  dens_id_entr = 1;
  active_id_entr = 1;
  dens_name[dens_id_entr] = "S";
  dens_desc[dens_id_entr] = "entropic variable density";
  dens_prognostic(dens_id_entr) = true;
  dens_active(dens_id_entr) = true;
  dens_pos(dens_id_entr) = false;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}

template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_entropic_var(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  return densvar(dens_id_entr, k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_entropic_var(const real3d &densvar,
                                                          int k, int ks,
                                                          int n) const {
  return densvar(dens_id_entr, k + ks, n) / densvar(dens_id_mass, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_alpha(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {
  return dual_geometry.get_area_n1entity(k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_CE>::get_alpha(const real3d &densvar, int k,
                                                   int ks, int n) const {
  return dual_geometry.get_area_n1entity(k + ks, 0, 0, n) /
         densvar(dens_id_mass, k + ks, n);
}
#endif

#ifdef _AN
template <>
void VariableSetBase<VS_AN>::initialize(PamCoupler &coupler,
                                        ModelParameters &params,
                                        const ThermoPotential &thermodynamics,
                                        const ReferenceState &refstate,
                                        const Geometry<Straight> &primal_geom,
                                        const Geometry<Twisted> &dual_geom,
                                        bool verbose) {
  dens_id_entr = 0;
  active_id_entr = 0;
  dens_name[dens_id_entr] = "S";
  dens_desc[dens_id_entr] = "entropic variable density";
  dens_prognostic(dens_id_entr) = true;
  dens_active(dens_id_entr) = true;
  dens_pos(dens_id_entr) = false;

  //  mass density is always stored last for anelastic
  //  this is to make prognostic densities a continuous subset of all densities
  dens_id_mass = ndensity - 1;
  active_id_mass = ndensity_active - 1;
  dens_name[dens_id_mass] = "rho";
  dens_desc[dens_id_mass] = "fluid density";
  dens_prognostic(dens_id_mass) = false;
  dens_active(dens_id_mass) = true;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_total_density(
    const real3d &densvar, int k, int ks, int n) const {
  return densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_entropic_var(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  return densvar(dens_id_entr, k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_entropic_var(const real3d &densvar,
                                                          int k, int ks,
                                                          int n) const {
  return densvar(dens_id_entr, k + ks, n) / densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_alpha(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {
  return dual_geometry.get_area_n1entity(k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_alpha(const real3d &densvar, int k,
                                                   int ks, int n) const {
  return dual_geometry.get_area_n1entity(k + ks, 0, 0, n) /
         densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_pres(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {
    const real refrho = 1.0_fp / get_alpha(reference_state.dens.data, k, ks, n);
    const real refentropic_var = get_entropic_var(reference_state.dens.data, k, ks, n);
    const real refp =  thermo.solve_p(refrho, refentropic_var, 0, 0, 0, 0);
    return refp;
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_pres(int k, int ks, int n) const {
    const real refrho = 1.0_fp / get_alpha(reference_state.dens.data, k, ks, n);
    const real refentropic_var = get_entropic_var(reference_state.dens.data, k, ks, n);
    const real refp =  thermo.solve_p(refrho, refentropic_var, 0, 0 , 0, 0);
    return refp;
}

template <>
real YAKL_INLINE VariableSetBase<VS_AN>::get_ref_dens(int k, int ks, int n) const {
    return reference_state.dens.data(dens_id_mass, k + ks, n);
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
                                         const Geometry<Twisted> &dual_geom,
                                         bool verbose) {
  dens_id_entr = 0;
  active_id_entr = 0;
  dens_name[dens_id_entr] = "S";
  dens_desc[dens_id_entr] = "entropic variable density";
  dens_prognostic(dens_id_entr) = true;
  dens_active(dens_id_entr) = true;
  dens_pos(dens_id_entr) = false;

  //  mass density is always stored last for anelastic
  //  this is to make prognostic densities a continuous subset of all densities
  dens_id_mass = ndensity - 1;
  active_id_mass = ndensity_active - 1;
  dens_name[dens_id_mass] = "rho";
  dens_desc[dens_id_mass] = "fluid density";
  dens_prognostic(dens_id_mass) = false;
  dens_active(dens_id_mass) = true;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_total_density(
    const real3d &densvar, int k, int ks, int n) const {
  return reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_entropic_var(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_entr, k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_entropic_var(
    const real3d &densvar, int k, int ks, int n) const {
  return reference_state.dens.data(dens_id_entr, k + ks, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_alpha(const real5d &densvar,
                                                    int k, int j, int i, int ks,
                                                    int js, int is,
                                                    int n) const {
  return dual_geometry.get_area_n1entity(k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_alpha(const real3d &densvar,
                                                    int k, int ks,
                                                    int n) const {
  return dual_geometry.get_area_n1entity(k + ks, 0, 0, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_ref_dens(int k, int ks, int n) const {
    return reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qv(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return densvar(dens_id_vap, k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qv(const real3d &densvar, int k,
                                                 int ks, int n) const {
  return densvar(dens_id_vap, k + ks, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_ql(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return densvar(dens_id_liq, k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qi(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return densvar(dens_id_ice, k + ks, j + js, i + is, n) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_pres(const real5d &densvar, int k,
                                                   int j, int i, int ks, int js,
                                                   int is, int n) const {
    const real refrho = 1.0_fp / get_alpha(reference_state.dens.data, k, ks, n);
    const real refentropic_var = get_entropic_var(reference_state.dens.data, k, ks, n);
    const real refqv = get_qv(reference_state.dens.data, k, ks, n);
    const real refp =  thermo.solve_p(refrho, refentropic_var, 1 - refqv, refqv, 0, 0);
    return refp;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_pres(int k, int ks, int n) const {
    const real refrho = 1.0_fp / get_alpha(reference_state.dens.data, k, ks, n);
    const real refentropic_var = get_entropic_var(reference_state.dens.data, k, ks, n);
    const real refqv = get_qv(reference_state.dens.data, k, ks, n);
    const real refp =  thermo.solve_p(refrho, refentropic_var, 1 - refqv, refqv, 0, 0);
    return refp;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::_water_dens(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  real vap_dens = densvar(dens_id_vap, k + ks, j + js, i + is, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  // if (liq_found) { liq_dens = densvar(dens_id_liq, k + ks, j + js, i + is, n); }
  // if (ice_found) { ice_dens = densvar(dens_id_ice, k + ks, j + js, i + is, n); }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::_water_dens(const real3d &densvar,
                                                      int k, int ks,
                                                      int n) const {
  real vap_dens = densvar(dens_id_vap, k + ks, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  // if (liq_found) { liq_dens = densvar(dens_id_liq, k + ks, n); }
  // if (ice_found) { ice_dens = densvar(dens_id_ice, k + ks, n); }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_dry_density(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  return (reference_state.dens.data(dens_id_mass, k + ks, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n));
}
template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qd(const real5d &densvar, int k,
                                                 int j, int i, int ks, int js,
                                                 int is, int n) const {
  return (reference_state.dens.data(dens_id_mass, k + ks, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n)) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MAN>::get_qd(const real3d &densvar, int k,
                                                 int ks, int n) const {
  return (reference_state.dens.data(dens_id_mass, k + ks, n) -
          _water_dens(densvar, k, ks, n)) /
         reference_state.dens.data(dens_id_mass, k + ks, n);
}

template <>
void YAKL_INLINE VariableSetBase<VS_MAN>::set_density(real dens, real dryden,
                                                      const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  reference_state.dens.data(dens_id_mass, k + ks, n) = dens;
}
template <>
void YAKL_INLINE VariableSetBase<VS_MAN>::set_entropic_density(
    real entropic_var_density, const real5d &densvar, int k, int j, int i,
    int ks, int js, int is, int n) const {
  densvar(dens_id_entr, k + ks, j + js, i + is, n) = entropic_var_density;
}
#endif

#ifdef _MCErho
template <>
void VariableSetBase<VS_MCE_rho>::initialize(
    PamCoupler &coupler, ModelParameters &params,
    const ThermoPotential &thermodynamics, const ReferenceState &refstate,
    const Geometry<Straight> &primal_geom, const Geometry<Twisted> &dual_geom,
    bool verbose) {
  dens_id_mass = 0;
  active_id_mass = 0;
  dens_name[dens_id_mass] = "rho";
  dens_desc[dens_id_mass] = "fluid density";
  dens_prognostic(dens_id_mass) = true;
  dens_active(dens_id_mass) = true;
  dens_pos(dens_id_mass) = false;

  dens_id_entr = 1;
  active_id_entr = 1;
  dens_desc[dens_id_entr] = "entropic variable density";
  dens_name[dens_id_entr] = "S";
  dens_prognostic(dens_id_entr) = true;
  dens_active(dens_id_entr) = true;
  dens_pos(dens_id_entr) = false;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_total_density(
    const real3d &densvar, int k, int ks, int n) const {
  return densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_entropic_var(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_entr, k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_entropic_var(
    const real3d &densvar, int k, int ks, int n) const {
  return densvar(dens_id_entr, k + ks, n) / densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_alpha(const real5d &densvar,
                                                        int k, int j, int i,
                                                        int ks, int js, int is,
                                                        int n) const {
  return dual_geometry.get_area_n1entity(k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_alpha(const real3d &densvar,
                                                        int k, int ks,
                                                        int n) const {
  return dual_geometry.get_area_n1entity(k + ks, 0, 0, n) /
         densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qv(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return densvar(dens_id_vap, k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qv(const real3d &densvar,
                                                     int k, int ks,
                                                     int n) const {
  return densvar(dens_id_vap, k + ks, n) / densvar(dens_id_mass, k + ks, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_ql(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return densvar(dens_id_liq, k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qi(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return densvar(dens_id_ice, k + ks, j + js, i + is, n) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::_water_dens(const real5d &densvar,
                                                          int k, int j, int i,
                                                          int ks, int js,
                                                          int is, int n) const {
  real vap_dens = densvar(dens_id_vap, k + ks, j + js, i + is, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liq_found) {
    liq_dens = densvar(dens_id_liq, k + ks, j + js, i + is, n);
  }
  if (ice_found) {
    ice_dens = densvar(dens_id_ice, k + ks, j + js, i + is, n);
  }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::_water_dens(const real3d &densvar,
                                                          int k, int ks,
                                                          int n) const {
  real vap_dens = densvar(dens_id_vap, k + ks, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liq_found) {
    liq_dens = densvar(dens_id_liq, k + ks, n);
  }
  if (ice_found) {
    ice_dens = densvar(dens_id_ice, k + ks, n);
  }
  return vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_dry_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return (densvar(dens_id_mass, k + ks, j + js, i + is, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n));
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qd(const real5d &densvar,
                                                     int k, int j, int i,
                                                     int ks, int js, int is,
                                                     int n) const {
  return (densvar(dens_id_mass, k + ks, j + js, i + is, n) -
          _water_dens(densvar, k, j, i, ks, js, is, n)) /
         densvar(dens_id_mass, k + ks, j + js, i + is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rho>::get_qd(const real3d &densvar,
                                                     int k, int ks,
                                                     int n) const {
  return (densvar(dens_id_mass, k + ks, n) - _water_dens(densvar, k, ks, n)) /
         densvar(dens_id_mass, k + ks, n);
}

template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rho>::set_density(
    real dens, real dryden, const real5d &densvar, int k, int j, int i, int ks,
    int js, int is, int n) const {
  densvar(dens_id_mass, k + ks, j + js, i + is, n) = dens;
}
template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rho>::set_entropic_density(
    real entropic_var_density, const real5d &densvar, int k, int j, int i,
    int ks, int js, int is, int n) const {
  densvar(dens_id_entr, k + ks, j + js, i + is, n) = entropic_var_density;
}
#endif

#ifdef _MCErhod
template <>
void VariableSetBase<VS_MCE_rhod>::initialize(
    PamCoupler &coupler, ModelParameters &params,
    const ThermoPotential &thermodynamics, const ReferenceState &refstate,
    const Geometry<Straight> &primal_geom, const Geometry<Twisted> &dual_geom,
    bool verbose) {
  dens_id_mass = 0;
  active_id_mass = 0;
  dens_name[dens_id_mass] = "rho_d";
  dens_desc[dens_id_mass] = "fluid dry density";
  dens_prognostic(dens_id_mass) = true;
  dens_active(dens_id_mass) = true;
  dens_pos(dens_id_mass) = false;

  dens_id_entr = 1;
  active_id_entr = 1;
  dens_name[dens_id_entr] = "S";
  dens_desc[dens_id_entr] = "entropic variable density";
  dens_prognostic(dens_id_entr) = true;
  dens_active(dens_id_entr) = true;
  dens_pos(dens_id_entr) = false;

  VariableSetBase::initialize(*this, coupler, params, thermodynamics, refstate,
                              primal_geom, dual_geom, verbose);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_total_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  real dry_dens = densvar(dens_id_mass, k + ks, j + js, i + is, n);
  real vap_dens = densvar(dens_id_vap, k + ks, j + js, i + is, n);
  real liq_dens = 0.0_fp;
  real ice_dens = 0.0_fp;
  if (liq_found) {
    liq_dens = densvar(dens_id_liq, k + ks, j + js, i + is, n);
  }
  if (ice_found) {
    ice_dens = densvar(dens_id_ice, k + ks, j + js, i + is, n);
  }
  return dry_dens + vap_dens + liq_dens + ice_dens;
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_dry_density(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_mass, k + ks, j + js, i + is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_entropic_var(
    const real5d &densvar, int k, int j, int i, int ks, int js, int is,
    int n) const {
  return densvar(dens_id_entr, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_alpha(const real5d &densvar,
                                                         int k, int j, int i,
                                                         int ks, int js, int is,
                                                         int n) const {
  return dual_geometry.get_area_n1entity(k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_qd(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dens_id_mass, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}

template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_qv(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dens_id_vap, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_ql(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dens_id_liq, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
real YAKL_INLINE VariableSetBase<VS_MCE_rhod>::get_qi(const real5d &densvar,
                                                      int k, int j, int i,
                                                      int ks, int js, int is,
                                                      int n) const {
  return densvar(dens_id_ice, k + ks, j + js, i + is, n) /
         get_total_density(densvar, k, j, i, ks, js, is, n);
}
template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rhod>::set_density(
    real dens, real drydens, const real5d &densvar, int k, int j, int i, int ks,
    int js, int is, int n) const {
  densvar(dens_id_mass, k + ks, j + js, i + is, n) = drydens;
}
template <>
void YAKL_INLINE VariableSetBase<VS_MCE_rhod>::set_entropic_density(
    real entropic_var_density, const real5d &densvar, int k, int j, int i,
    int ks, int js, int is, int n) const {
  densvar(dens_id_entr, k + ks, j + js, i + is, n) = entropic_var_density;
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
 
