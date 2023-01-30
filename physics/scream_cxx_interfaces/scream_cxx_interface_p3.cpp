
#include "EKAT_utils.h"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"
#include <fstream>      // std::ifstream

void p3_main_cxx(array_ir::ArrayIR<double,2> const & qc,                 // inout
                 array_ir::ArrayIR<double,2> const & nc,                 // inout
                 array_ir::ArrayIR<double,2> const & qr,                 // inout
                 array_ir::ArrayIR<double,2> const & nr,                 // inout
                 array_ir::ArrayIR<double,2> const & th_atm,             // inout
                 array_ir::ArrayIR<double,2> const & qv,                 // inout
                 double                      const & dt,                 // inout
                 array_ir::ArrayIR<double,2> const & qi,                 // inout
                 array_ir::ArrayIR<double,2> const & qm,                 // inout
                 array_ir::ArrayIR<double,2> const & ni,                 // inout
                 array_ir::ArrayIR<double,2> const & bm,                 // inout
                 array_ir::ArrayIR<double,2> const & pres,               // in
                 array_ir::ArrayIR<double,2> const & dz,                 // in
                 array_ir::ArrayIR<double,2> const & nc_nuceat_tend,     // in
                 array_ir::ArrayIR<double,2> const & nccn_prescribed,    // in
                 array_ir::ArrayIR<double,2> const & ni_activated,       // in
                 array_ir::ArrayIR<double,2> const & inv_qc_relvar,      // in
                 int                         const & it,                 // in
                 array_ir::ArrayIR<double,1> const & precip_liq_surf,    //   out
                 array_ir::ArrayIR<double,1> const & precip_ice_surf,    //   out
                 int                         const & its,                // in
                 int                         const & ite,                // in
                 int                         const & kts,                // in
                 int                         const & kte,                // in
                 array_ir::ArrayIR<double,2> const & diag_eff_radius_qc, //   out
                 array_ir::ArrayIR<double,2> const & diag_eff_radius_qi, //   out
                 array_ir::ArrayIR<double,2> const & rho_qi,             //   out
                 bool                        const & do_predict_nc,      // in
                 bool                        const & do_prescribed_CCN,  // in
                 array_ir::ArrayIR<double,2> const & dpres,              // in
                 array_ir::ArrayIR<double,2> const & inv_exner,          // in
                 array_ir::ArrayIR<double,2> const & qv2qi_depos_tend,   //   out
                 array_ir::ArrayIR<double,2> const & precip_liq_flux,    //   out
                 array_ir::ArrayIR<double,2> const & precip_ice_flux,    //   out
                 array_ir::ArrayIR<double,2> const & cld_frac_r,         // in
                 array_ir::ArrayIR<double,2> const & cld_frac_l,         // in
                 array_ir::ArrayIR<double,2> const & cld_frac_i,         // in
                 array_ir::ArrayIR<double,2> const & liq_ice_exchange,   //   out
                 array_ir::ArrayIR<double,2> const & vap_liq_exchange,   //   out
                 array_ir::ArrayIR<double,2> const & vap_ice_exchange,   //   out
                 array_ir::ArrayIR<double,2> const & qv_prev,            // in
                 array_ir::ArrayIR<double,2> const & t_prev,             // in
                 array_ir::ArrayIR<double,2> const & col_location,       // in
                 double                              *elapsed_s ) {      //   out
  using ScreamCXX::ArrayIR_to_View_of_Packs;
  using ScreamCXX::ArrayIR_to_View;
  using namespace scream;
  using namespace scream::p3;
  using P3F        = p3::Functions<Real, DefaultDevice>;
  using KT         = typename P3F::KT;
  using Spack      = typename P3F::Spack;

  const int nlev  = kte-kts+1;
  const int ncol  = ite-its+1;
  const int npack = ekat::npack<Spack>(nlev);

  auto qv_d = ArrayIR_to_View_of_Packs(qv);
  auto qc_d = ArrayIR_to_View_of_Packs(qc);
  auto nc_d = ArrayIR_to_View_of_Packs(nc);
  auto qr_d = ArrayIR_to_View_of_Packs(qr);
  auto nr_d = ArrayIR_to_View_of_Packs(nr);
  auto qi_d = ArrayIR_to_View_of_Packs(qi);
  auto qm_d = ArrayIR_to_View_of_Packs(qm);
  auto ni_d = ArrayIR_to_View_of_Packs(ni);
  auto bm_d = ArrayIR_to_View_of_Packs(bm);
  auto th_d = ArrayIR_to_View_of_Packs(th_atm);

  P3F::P3PrognosticState prog_state{qc_d, nc_d, qr_d, nr_d, qi_d, qm_d,
                                    ni_d, bm_d, qv_d, th_d};

  auto nc_nuceat_tend_d = ArrayIR_to_View_of_Packs(nc_nuceat_tend );
  auto nccn_d           = ArrayIR_to_View_of_Packs(nccn_prescribed);
  auto ni_activated_d   = ArrayIR_to_View_of_Packs(ni_activated   );
  auto inv_qc_relvar_d  = ArrayIR_to_View_of_Packs(inv_qc_relvar  );
  auto dz_d             = ArrayIR_to_View_of_Packs(dz             );
  auto pres_d           = ArrayIR_to_View_of_Packs(pres           );
  auto dpres_d          = ArrayIR_to_View_of_Packs(dpres          );
  auto inv_exner_d      = ArrayIR_to_View_of_Packs(inv_exner      );
  auto t_prev_d         = ArrayIR_to_View_of_Packs(t_prev         );
  auto q_prev_d         = ArrayIR_to_View_of_Packs(qv_prev        );
  auto cld_frac_i_d     = ArrayIR_to_View_of_Packs(cld_frac_i     );
  auto cld_frac_l_d     = ArrayIR_to_View_of_Packs(cld_frac_l     );
  auto cld_frac_r_d     = ArrayIR_to_View_of_Packs(cld_frac_r     );

  P3F::P3DiagnosticInputs diag_inputs{nc_nuceat_tend_d, nccn_d, 
                                      ni_activated_d, inv_qc_relvar_d, 
                                      cld_frac_i_d, cld_frac_l_d, cld_frac_r_d, 
                                      pres_d, dz_d, dpres_d, inv_exner_d, 
                                      q_prev_d, t_prev_d};

  auto qv2qi_depos_tend_d   = ArrayIR_to_View_of_Packs(qv2qi_depos_tend  );
  auto diag_eff_radius_qc_d = ArrayIR_to_View_of_Packs(diag_eff_radius_qc);
  auto diag_eff_radius_qi_d = ArrayIR_to_View_of_Packs(diag_eff_radius_qi);
  auto rho_qi_d             = ArrayIR_to_View_of_Packs(rho_qi            );
  auto precip_liq_flux_d    = ArrayIR_to_View_of_Packs(precip_liq_flux   );
  auto precip_ice_flux_d    = ArrayIR_to_View_of_Packs(precip_ice_flux   );
  auto precip_liq_surf_d    = ArrayIR_to_View         (precip_liq_surf   );
  auto precip_ice_surf_d    = ArrayIR_to_View         (precip_ice_surf   );

  Kokkos::parallel_for( Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}),
                        KOKKOS_LAMBDA(int icol, int ilev, int s) {
    qv2qi_depos_tend_d  (icol,ilev)[s] = 0.;
    diag_eff_radius_qc_d(icol,ilev)[s] = 0.;
    diag_eff_radius_qi_d(icol,ilev)[s] = 0.;
    rho_qi_d            (icol,ilev)[s] = 0.;
    precip_liq_flux_d   (icol,ilev)[s] = 0.;
    precip_ice_flux_d   (icol,ilev)[s] = 0.;
    if (ilev == 0 && s == 0) {
      precip_liq_surf_d(icol) = 0.;
      precip_ice_surf_d(icol) = 0.;
    }
  });

  P3F::P3DiagnosticOutputs diag_outputs {qv2qi_depos_tend_d, precip_liq_surf_d,
                                         precip_ice_surf_d, diag_eff_radius_qc_d, diag_eff_radius_qi_d,
                                         rho_qi_d,precip_liq_flux_d, precip_ice_flux_d};

  auto col_location_d = ArrayIR_to_View(col_location);

  P3F::P3Infrastructure infrastructure{dt, it, 0, ite-its, 0, kte-kts,
                                       do_predict_nc, do_prescribed_CCN, col_location_d};

  auto liq_ice_exchange_d = ArrayIR_to_View_of_Packs(liq_ice_exchange);
  auto vap_liq_exchange_d = ArrayIR_to_View_of_Packs(vap_liq_exchange);
  auto vap_ice_exchange_d = ArrayIR_to_View_of_Packs(vap_ice_exchange);

  Kokkos::parallel_for( Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, npack}),
                        KOKKOS_LAMBDA(int icol, int ilev) {
     liq_ice_exchange_d(icol,ilev) = 0.;
     vap_liq_exchange_d(icol,ilev) = 0.;
     vap_ice_exchange_d(icol,ilev) = 0.;
  });

  P3F::P3HistoryOnly history_only {liq_ice_exchange_d, vap_liq_exchange_d, vap_ice_exchange_d};

  P3F::view_1d_table      mu_r_table_vals;
  P3F::view_2d_table      vn_table_vals, vm_table_vals, revap_table_vals;
  P3F::view_ice_table     ice_table_vals;
  P3F::view_collect_table collect_table_vals;
  P3F::view_dnu_table     dnu_table_vals;
  P3F::init_kokkos_ice_lookup_tables(ice_table_vals, collect_table_vals);
  P3F::init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu_table_vals);

  P3F::P3LookupTables tables{mu_r_table_vals, vn_table_vals, vm_table_vals, revap_table_vals,
                             ice_table_vals, collect_table_vals, dnu_table_vals};

  const int nlev_pack = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol, nlev_pack);
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(nlev_pack, 59, policy);

  auto elapsed_time = P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
                                   history_only, tables, workspace_mgr, ncol, nlev);
}


