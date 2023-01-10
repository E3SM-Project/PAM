
#pragma once

#include "pam_utils.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

#include <fstream>      // std::ifstream

#define USE_SCREAM

using namespace scream;
using namespace scream::p3;

using array1D = IRType<Real, 1>::type;
using array2D = IRType<Real, 2>::type;
using array3D = IRType<Real, 3>::type;

//
// main micro_p3 microproc
//
void p3_main_cxx(array2D& qc, array2D& nc, array2D& qr, array2D& nr, array2D& th_atm, array2D&qv,
                 double &dt, array2D& qi, array2D& qm, array2D& ni, array2D& bm, array2D& pres,
                 array2D& dz, array2D& nc_nuceat_tend, array2D& nccn_prescribed, array2D& ni_activated,
                 array2D& inv_qc_relvar, int &it, array1D& precip_liq_surf, array1D& precip_ice_surf,
                 int &its, int &ite, int &kts, int &kte, array2D& diag_eff_radius_qc,
                 array2D& diag_eff_radius_qi, array2D& rho_qi, bool &do_predict_nc,
                 bool &do_prescribed_CCN, array2D& dpres, array2D& exner, array2D& qv2qi_depos_tend,
                 array2D& precip_total_tend, array2D& nevapr, array2D& qr_evap_tend,
                 array2D& precip_liq_flux, array2D& precip_ice_flux, array2D& cld_frac_r,
                 array2D& cld_frac_l, array2D& cld_frac_i, array3D& p3_tend_out, array2D& mu_c,
                 array2D& lamc, array2D& liq_ice_exchange, array2D& vap_liq_exchange,
                 array2D& vap_ice_exchange, array2D& qv_prev, array2D& t_prev, array2D& col_location,
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

  auto qv_in = reshape(qv);
  auto qc_in = reshape(qc);
  auto nc_in = reshape(nc);
  auto qr_in = reshape(qr);
  auto nr_in = reshape(nr);
  auto qi_in = reshape(qi);
  auto qm_in = reshape(qm);
  auto ni_in = reshape(ni);
  auto bm_in = reshape(bm);
  auto th_in = reshape(th_atm);

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
  auto nc_nuceat_tend_in = reshape(nc_nuceat_tend);
  auto nccn_in           = reshape(nccn_prescribed);
  auto ni_activated_in   = reshape(ni_activated);
  auto inv_qc_relvar_in  = reshape(inv_qc_relvar);
  auto dz_in             = reshape(dz);
  auto pres_in           = reshape(pres);
  auto dpres_in          = reshape(dpres);
  auto t_prev_in         = reshape(t_prev);
  auto q_prev_in         = reshape(qv_prev);
  auto cld_frac_i_in     = reshape(cld_frac_i);
  auto cld_frac_l_in     = reshape(cld_frac_l);
  auto cld_frac_r_in     = reshape(cld_frac_r);

  auto inv_exner_in      = reshape_and_inverse(exner);

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
     col_location_d(i, 0) = col_location.data()[i*3+0];
     col_location_d(i, 1) = col_location.data()[i*3+1];
     col_location_d(i, 2) = col_location.data()[i*3+2];
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
    int k = ilev*Spack::n + s;
    int offset = k*ncol + icol;
    if (k < nlev) {
      qv.data()[offset] = prog_state.qv(icol,ilev)[s] + prog_state.qc(icol,ilev)[s];
      qc.data()[offset] = prog_state.qc(icol,ilev)[s];
      nc.data()[offset] = prog_state.nc(icol,ilev)[s];
      qr.data()[offset] = prog_state.qr(icol,ilev)[s];
      nr.data()[offset] = prog_state.nr(icol,ilev)[s];
      qi.data()[offset] = prog_state.qi(icol,ilev)[s];
      qm.data()[offset] = prog_state.qm(icol,ilev)[s];
      ni.data()[offset] = prog_state.ni(icol,ilev)[s];
      bm.data()[offset] = prog_state.bm(icol,ilev)[s];
    }
  });

  // update LSE, temperature, and previous t/q
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int k    = ilev*Spack::n + s;
    int offset = k*ncol + icol;
    if (k < nlev) {
      t_prev.data()[offset]  = prog_state.th(icol,ilev)[s]/inv_exner_d(icol,ilev)[s];
      qv_prev.data()[offset] = qv.data()[icol, ilev];
    }
  });
}

