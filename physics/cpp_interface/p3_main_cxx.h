
#pragma once

#include "pam_utils.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

#include <fstream>      // std::ifstream

#define USE_SCREAM

using namespace scream;
using namespace scream::p3;

// reshape an input array (nj, nk) to output array (nk, nj)
template <typename Scalar>
void reshape(const Scalar* sv, Scalar* dv, Int nj, Int nk) {
  parallel_for( SimpleBounds<2>(nj, nk) , YAKL_LAMBDA (int j, int k) {
    dv[k*nj+j] = sv[j*nk+k];
  });
}

// reshape and inverse an input array (nj, nk) to output array (nk, nj)
template <typename Scalar>
void reshape_and_inverse(const Scalar* sv, Scalar* dv, Int nj, Int nk) {
  parallel_for( SimpleBounds<2>(nj, nk) , YAKL_LAMBDA (int j, int k) {
    dv[k*nj+j] = 1./sv[j*nk+k];
  });
}

//
// main micro_p3 microproc
//
void p3_main_cxx(double *qc , double *nc , double *qr , double *nr , double *th_atm , double *qv ,
                 double &dt , double *qi , double *qm , double *ni , double *bm , double *pres ,
                 double *dz , double *nc_nuceat_tend , double *nccn_prescribed , double *ni_activated ,
                 double *inv_qc_relvar , int &it , double *precip_liq_surf , double *precip_ice_surf ,
                 int &its , int &ite , int &kts , int &kte , double *diag_eff_radius_qc ,
                 double *diag_eff_radius_qi , double *rho_qi , bool &do_predict_nc ,
                 bool &do_prescribed_CCN ,double *dpres , double *exner , double *qv2qi_depos_tend ,
                 double *precip_total_tend , double *nevapr , double *qr_evap_tend ,
                 double *precip_liq_flux , double *precip_ice_flux , double *cld_frac_r ,
                 double *cld_frac_l , double *cld_frac_i , double *p3_tend_out , double *mu_c ,
                 double *lamc , double *liq_ice_exchange , double *vap_liq_exchange ,
                 double *vap_ice_exchange , double *qv_prev , double *t_prev , double *col_location ,
                 double *elapsed_s )
{
  using P3F        = p3::Functions<Real, DefaultDevice>;
  using KT         = typename P3F::KT;
  using ExeSpace   = typename KT::ExeSpace;
  using MemberType = typename P3F::MemberType;
  using Spack      = typename P3F::Spack;
  using view_1d    = typename P3F::view_1d<Spack>;
  using uview_1d   = typename P3F::uview_1d<Spack>;
  using view_2d    = typename P3F::view_2d<Spack>;
  using sview_1d   = typename P3F::view_1d<Real>;
  using sview_2d   = typename P3F::view_2d<Real>;

  using view_1d_table      = typename P3F::view_1d_table;
  using view_2d_table      = typename P3F::view_2d_table;
  using view_ice_table     = typename P3F::view_ice_table;
  using view_collect_table = typename P3F::view_collect_table;
  using view_dnu_table     = typename P3F::view_dnu_table;

  const int nlev  = kte-kts+1;
  const int ncol  = ite-its+1;
  const int npack = ekat::npack<Spack>(nlev);

  real2d qc_in("qc",ncol, nlev);
  real2d nc_in("nc",ncol, nlev);
  real2d qr_in("qr",ncol, nlev);
  real2d nr_in("nr",ncol, nlev);
  real2d qi_in("qi",ncol, nlev);
  real2d qm_in("qm",ncol, nlev);
  real2d ni_in("ni",ncol, nlev);
  real2d bm_in("bm",ncol, nlev);
  real2d qv_in("qv",ncol, nlev);
  real2d th_in("th",ncol, nlev);

  real2d nc_nuceat_tend_in("nuceat",ncol, nlev);
  real2d nccn_in("nccn",ncol, nlev);
  real2d ni_activated_in("ni_act",ncol, nlev);
  real2d inv_qc_relvar_in("inv_qc",ncol, nlev); // relative cloud water variance - not needed - set to 1
  real2d cld_frac_i_in("cld_frac_i",ncol, nlev);
  real2d cld_frac_l_in("cld_frac_l",ncol, nlev);
  real2d cld_frac_r_in("cld_frac_r",ncol, nlev);
  real2d dz_in("dz", ncol, nlev);
  real2d pres_in("pres", ncol, nlev);
  real2d dpres_in("dpres",ncol, nlev);
  real2d inv_exner_in("inv_exner",ncol, nlev);
  real2d q_prev_in("q_prev",ncol, nlev);
  real2d t_prev_in("t_prev",ncol, nlev);

  real2d cloud_frac_in("cloud_frac_in",ncol, nlev);

  reshape(qv,     qv_in.data(), nlev, ncol);
  reshape(qc,     qc_in.data(), nlev, ncol);
  reshape(nc,     nc_in.data(), nlev, ncol);
  reshape(qr,     qr_in.data(), nlev, ncol);
  reshape(nr,     nr_in.data(), nlev, ncol);
  reshape(qi,     qi_in.data(), nlev, ncol);
  reshape(qm,     qm_in.data(), nlev, ncol);
  reshape(ni,     ni_in.data(), nlev, ncol);
  reshape(bm,     bm_in.data(), nlev, ncol);
  reshape(th_atm, th_in.data(), nlev, ncol);

  view_2d qv_d("qv", ncol, npack),
          qc_d("qc", ncol, npack),
          nc_d("nc", ncol, npack),
          qr_d("qr", ncol, npack),
          nr_d("nr", ncol, npack),
          qi_d("qi", ncol, npack),
          qm_d("qm", ncol, npack),
          ni_d("ni", ncol, npack),
          bm_d("bm", ncol, npack),
          th_d("th", ncol, npack);

  array_to_view(qc_in.data(), ncol, nlev, qc_d);
  array_to_view(nc_in.data(), ncol, nlev, nc_d);
  array_to_view(qr_in.data(), ncol, nlev, qr_d);
  array_to_view(nr_in.data(), ncol, nlev, nr_d);
  array_to_view(qi_in.data(), ncol, nlev, qi_d);
  array_to_view(qm_in.data(), ncol, nlev, qm_d);
  array_to_view(ni_in.data(), ncol, nlev, ni_d);
  array_to_view(bm_in.data(), ncol, nlev, bm_d);
  array_to_view(qv_in.data(), ncol, nlev, qv_d);
  array_to_view(th_in.data(), ncol, nlev, th_d);
   
  P3F::P3PrognosticState prog_state{qc_d, nc_d, qr_d, nr_d, qi_d, qm_d,
                                    ni_d, bm_d, qv_d, th_d};

  //----------------------------------------------------------------------------
  // Populate P3 diagnostic inputs
  //----------------------------------------------------------------------------
  reshape(nc_nuceat_tend, nc_nuceat_tend_in.data(), nlev, ncol);
  reshape(nccn_prescribed,nccn_in.data(),           nlev, ncol);
  reshape(ni_activated,   ni_activated_in.data(),   nlev, ncol);
  reshape(inv_qc_relvar,  inv_qc_relvar_in.data(),  nlev, ncol);
  reshape(dz,             dz_in.data(),             nlev, ncol);
  reshape(pres,           pres_in.data(),           nlev, ncol);
  reshape(dpres,          dpres_in.data(),          nlev, ncol);
  reshape(t_prev,         t_prev_in.data(),         nlev, ncol);
  reshape(qv_prev,        q_prev_in.data(),         nlev, ncol);
  reshape(cld_frac_i,     cld_frac_i_in.data(),     nlev, ncol);
  reshape(cld_frac_l,     cld_frac_l_in.data(),     nlev, ncol);
  reshape(cld_frac_r,     cld_frac_r_in.data(),     nlev, ncol);

  reshape_and_inverse(exner,inv_exner_in.data(), nlev, ncol);

  view_2d nc_nuceat_tend_d("nc_nuceat_tend", ncol, npack),
          nccn_d("nccn", ncol, npack),
          ni_activated_d("ni_activated", ncol, npack),
          inv_qc_relvar_d("inv_qc_relvar", ncol, npack),
          dz_d("dz", ncol, npack),
          pres_d("pres", ncol, npack),
          dpres_d("dpres", ncol, npack),
          inv_exner_d("inv_exner", ncol, npack),
          t_prev_d("t_prev", ncol, npack),
          q_prev_d("q_prev", ncol, npack),
          cld_frac_i_d("cld_frac_i", ncol, npack),
          cld_frac_l_d("cld_frac_l", ncol, npack),
          cld_frac_r_d("cld_frac_r", ncol, npack);

  array_to_view(nc_nuceat_tend_in.data(), ncol, nlev, nc_nuceat_tend_d);
  array_to_view(nccn_in.data(), ncol, nlev, nccn_d);
  array_to_view(ni_activated_in.data(), ncol, nlev, ni_activated_d);
  array_to_view(inv_qc_relvar_in.data(), ncol, nlev, inv_qc_relvar_d);
  array_to_view(dz_in.data(), ncol, nlev, dz_d);
  array_to_view(pres_in.data(), ncol, nlev, pres_d);
  array_to_view(dpres_in.data(), ncol, nlev, dpres_d);
  array_to_view(inv_exner_in.data(), ncol, nlev, inv_exner_d);
  array_to_view(t_prev_in.data(), ncol, nlev, t_prev_d);
  array_to_view(q_prev_in.data(), ncol, nlev, q_prev_d);
  array_to_view(cld_frac_i_in.data(), ncol, nlev, cld_frac_i_d);
  array_to_view(cld_frac_l_in.data(), ncol, nlev, cld_frac_l_d);
  array_to_view(cld_frac_r_in.data(), ncol, nlev, cld_frac_r_d);

  P3F::P3DiagnosticInputs diag_inputs{nc_nuceat_tend_d, nccn_d, 
                                      ni_activated_d, inv_qc_relvar_d, 
                                      cld_frac_i_d, cld_frac_l_d, cld_frac_r_d, 
                                      pres_d, dz_d, dpres_d, inv_exner_d, 
                                      q_prev_d, t_prev_d};

  //----------------------------------------------------------------------------
  // Populate P3 diagnostic outputs
  //----------------------------------------------------------------------------
  view_2d qv2qi_depos_tend_d("qv2qi_depos_tend", ncol, npack),
          diag_eff_radius_qc_d("diag_eff_radius_qc", ncol, npack),
          diag_eff_radius_qi_d("diag_eff_radius_qi", ncol, npack),
          rho_qi_d("rho_qi", ncol, npack),
          precip_liq_flux_d("precip_liq_flux", ncol, npack),
          precip_ice_flux_d("precip_ice_flux", ncol, npack);

  sview_1d precip_liq_surf_d("precip_liq_surf_d", ncol), 
           precip_ice_surf_d("precip_ice_surf_d", ncol);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
     qv2qi_depos_tend_d(icol,ilev)[s] = 0.;
     diag_eff_radius_qc_d(icol,ilev)[s] = 0.;
     diag_eff_radius_qi_d(icol,ilev)[s] = 0.;
     rho_qi_d(icol,ilev)[s] = 0.;
     precip_liq_flux_d(icol,ilev)[s] = 0.;
     precip_ice_flux_d(icol,ilev)[s] = 0.;
  });

  Kokkos::parallel_for("precip", ncol, KOKKOS_LAMBDA (const int& icol) {
    precip_liq_surf_d(icol) = 0.;
    precip_ice_surf_d(icol) = 0.;
  });

  P3F::P3DiagnosticOutputs diag_outputs {qv2qi_depos_tend_d, precip_liq_surf_d,
                                         precip_ice_surf_d, diag_eff_radius_qc_d, diag_eff_radius_qi_d,
                                         rho_qi_d,precip_liq_flux_d, precip_ice_flux_d};

  //----------------------------------------------------------------------------
  // Populate P3 infrastructure
  //----------------------------------------------------------------------------
  sview_2d col_location_d("col_location_d", ncol, 3);
  Kokkos::parallel_for("col_location", ncol, KOKKOS_LAMBDA (const int& i) {
     col_location_d(i, 0) = col_location[i*3+0];
     col_location_d(i, 1) = col_location[i*3+1];
     col_location_d(i, 2) = col_location[i*3+2];
  });

  P3F::P3Infrastructure infrastructure{dt, it, 0, ite-its, 0, kte-kts,
                                       do_predict_nc, do_prescribed_CCN, col_location_d};

  //----------------------------------------------------------------------------
  // Populate P3 history output
  //----------------------------------------------------------------------------
  view_2d liq_ice_exchange_d("liq_ice_exchange_d", ncol, npack),
          vap_liq_exchange_d("vap_liq_exchange_d", ncol, npack),
          vap_ice_exchange_d("vap_ice_exchange_d", ncol, npack);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, npack}), KOKKOS_LAMBDA(int icol, int ilev) {
     liq_ice_exchange_d(icol,ilev) = 0.;
     vap_liq_exchange_d(icol,ilev) = 0.;
     vap_ice_exchange_d(icol,ilev) = 0.;
  });

  P3F::P3HistoryOnly history_only {liq_ice_exchange_d, vap_liq_exchange_d,
                                   vap_ice_exchange_d};


  // load tables
  view_1d_table mu_r_table_vals;
  view_2d_table vn_table_vals, vm_table_vals, revap_table_vals;
  view_ice_table ice_table_vals;
  view_collect_table collect_table_vals;
  view_dnu_table dnu_table_vals;
  P3F::init_kokkos_ice_lookup_tables(ice_table_vals, collect_table_vals);
  P3F::init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu_table_vals);

  P3F::P3LookupTables tables{mu_r_table_vals, vn_table_vals, vm_table_vals, revap_table_vals,
                             ice_table_vals, collect_table_vals, dnu_table_vals};

  //----------------------------------------------------------------------------
  // Call p3_main
  //----------------------------------------------------------------------------
  const int nlev_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol, nlev_pack);
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlev_pack, 59, policy);

  auto elapsed_time = P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
                                   history_only, tables, workspace_mgr, ncol, nlev);

//printf("p3_main wall time: nj=%d, nk=%d, time=%13.6e\n", ite, kte, (float)elapsed_time*1.e-6);
   
  // update microfield
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int k    = ilev*Spack::n + s;
    if (k < nlev) {
      qv_in(icol, ilev) = prog_state.qv(icol,ilev)[s] + prog_state.qc(icol,ilev)[s];
      qc_in(icol, ilev) = prog_state.qc(icol,ilev)[s];
      nc_in(icol, ilev) = prog_state.nc(icol,ilev)[s];
      qr_in(icol, ilev) = prog_state.qr(icol,ilev)[s];
      nr_in(icol, ilev) = prog_state.nr(icol,ilev)[s];
      qi_in(icol, ilev) = prog_state.qi(icol,ilev)[s];
      qm_in(icol, ilev) = prog_state.qm(icol,ilev)[s];
      ni_in(icol, ilev) = prog_state.ni(icol,ilev)[s];
      bm_in(icol, ilev) = prog_state.bm(icol,ilev)[s];
    }
  });

  // update LSE, temperature, and previous t/q
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int k    = ilev*Spack::n + s;
    if (k < nlev) {
      t_prev_in(icol, ilev) = prog_state.th(icol,ilev)[s]/inv_exner_in(icol,ilev);
      q_prev_in(icol, ilev) = qv_in(icol, ilev);
    }
  });

  reshape(qv_in.data(),     qv,     ncol, nlev);
  reshape(qc_in.data(),     qc,     ncol, nlev);
  reshape(nc_in.data(),     nc,     ncol, nlev);
  reshape(qr_in.data(),     qr,     ncol, nlev);
  reshape(nr_in.data(),     nr,     ncol, nlev);
  reshape(qi_in.data(),     qi,     ncol, nlev);
  reshape(qm_in.data(),     qm,     ncol, nlev);
  reshape(ni_in.data(),     ni,     ncol, nlev);
  reshape(bm_in.data(),     bm,     ncol, nlev);
  reshape(t_prev_in.data(), t_prev, ncol, nlev);
  reshape(q_prev_in.data(), qv,     ncol, nlev);

}

