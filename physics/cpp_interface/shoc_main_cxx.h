
#pragma once

#include "pam_utils.h"

#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "shoc_f90.hpp"

using namespace scream;
using namespace scream::shoc;

void shoc_main_cxx(int &shcol, int &nlev, int &nlevi, double &dtime, int &nadv, double *host_dx, double *host_dy,
                   double *thv, double *zt_grid, double *zi_grid, double *pres, double *presi, double *pdel,
                   double *wthl_sfc, double *wqw_sfc, double *uw_sfc, double *vw_sfc, double *wtracer_sfc,
                   int &num_qtracers, double *w_field, double *inv_exner, double *phis, double *host_dse, double *tke,
                   double *thetal, double *qw, double *u_wind, double *v_wind, double *qtracers, double *wthv_sec,
                   double *tkh, double *tk, double *shoc_ql, double *shoc_cldfrac, double *pblh, double *shoc_mix,
                   double *isotropy, double *w_sec, double *thl_sec, double *qw_sec, double *qwthl_sec,
                   double *wthl_sec, double *wqw_sec, double *wtke_sec, double *uw_sec, double *vw_sec, double *w3,
                   double *wqls_sec, double *brunt, double *shoc_ql2 )

{
  using SHOC       = shoc::Functions<Real, DefaultDevice>;
  using Scalar     = typename SHOC::Scalar;
  using Spack      = typename SHOC::Spack;
  using view_1d    = typename SHOC::view_1d<Scalar>;
  using view_2d    = typename SHOC::view_2d<Spack>;
  using view_3d    = typename SHOC::view_3d<Spack>;
  using ExeSpace   = typename SHOC::KT::ExeSpace;
  using MemberType = typename SHOC::MemberType;

  const int ncol   = shcol;
  const int npack  = ekat::npack<Spack>(nlev);
  const int nipack = ekat::npack<Spack>(nlevi);

  real2d inv_exner_in("inv_exner", ncol, nlev);
  real2d qw_in("qw", ncol, nlev); // total water (vapor + cloud liquid)
  real2d shoc_ql_in("shoc_ql", ncol, nlev); // cloud liquid water
  real2d thv_in("thv", ncol, nlev);
  real2d zt_g_in("zt_g", ncol, nlev);
  real2d zi_g_in("zi_g", ncol, nlevi);
  real2d tke_in("tke", ncol, nlev);
  real2d tk_in("tk", ncol, nlev);
  real2d tkh_in("tkh", ncol, nlev);
  real2d wthv_sec_in("wthv_sec", ncol, nlev);
  real2d pres_in("pres", ncol, nlev);
  real2d presi_in("presi", ncol, nlev);
  real2d pdel_in("pdel", ncol, nlevi);
  real2d host_dse_in("host_dse", ncol, nlev);
  real2d shoc_dse_in("shoc_dse", ncol, nlev);
  real2d shoc_cldfrac_in("shoc_cldfrac", ncol, nlev);
  real2d thetal_in("thetal", ncol, nlev);
  real2d w_field_in("w_field", ncol, nlev);
  real2d wtracer_sfc_in("wtracer_sfc", ncol, nlev);
  real3d shoc_hwind_in("shoc_hwind", ncol, 2, nlev);
  real3d qtracers_in("qtracer", ncol, num_qtracers, nlev);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}), KOKKOS_LAMBDA(int i, int k) {
     const int offset = i*nlev+k;
     shoc_hwind_in(i,0,k) = u_wind[offset];
     shoc_hwind_in(i,1,k) = v_wind[offset];
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, num_qtracers, nlev}), KOKKOS_LAMBDA(int i, int q, int k) {
    const int offset = (i*nlev+k)*num_qtracers+q;
    qtracers_in(i,q,k) = qtracers[offset];
  });

  // -------------------------------------------------
  // Set surface geopotential and fluxes
  // -------------------------------------------------
  view_1d host_dx_1d(host_dx),
          host_dy_1d(host_dy),
          wthl_sfc_1d(wthl_sfc),  // Surface sensible heat flux [K m/s]
          wqw_sfc_1d(wqw_sfc),    // Surface latent heat flux [kg/kg m/s]
          uw_sfc_1d(uw_sfc),      // Surface momentum flux (u-direction) [m2/s2]
          vw_sfc_1d(vw_sfc),      // Surface momentum flux (v-direction) [m2/s2]
          phis_1d(phis_1d);       // Host model surface geopotential height

  reshape(zt_grid,     zt_g_in.data(),       ncol, nlev);
  reshape(zi_grid,     zi_g_in.data(),       ncol, nlevi);
  reshape(pres,        pres_in.data(),       ncol, nlev);
  reshape(presi,       presi_in.data(),      ncol, nlevi);
  reshape(pdel,        pdel_in.data(),       ncol, nlev);
  reshape(inv_exner,   inv_exner_in.data(),  ncol, nlev);
  reshape(thv,         thv_in.data(),        ncol, nlev);
  reshape(w_field,     w_field_in.data(),    ncol, nlev);
  reshape(wtracer_sfc, wtracer_sfc_in.data(),ncol, nlev);

  view_2d zt_grid_2d("zt_grid", ncol, npack),                 // heights, for thermo grid [m]
          zi_grid_2d("zi_grid", ncol, nipack),                // heights, for interface grid [m]
          pres_2d("pres", ncol, npack),                       // pressure levels on thermo grid [Pa]
          presi_2d("presi", ncol, nipack),                    // pressure levels on interface grid [Pa]
          pdel_2d("pdel", ncol, npack),                       // Differences in pressure levels [Pa]
          thv_2d("thv", ncol, npack),                         // virtual potential temperature [K]
          w_field_2d("w_field", ncol, npack),                 // large scale vertical velocity [m/s]
          wtracer_sfc_2d("wtracer", ncol, npack),             // Surface flux for tracers [varies]
          inv_exner_2d("inv_exner", ncol, npack);

  array_to_view(zt_g_in.data(),      ncol, nlev,  zt_grid_2d);
  array_to_view(zi_g_in.data(),      ncol, nlevi, zi_grid_2d);
  array_to_view(pres_in.data(),      ncol, nlev,  pres_2d);
  array_to_view(presi_in.data(),     ncol, nlevi, presi_2d);
  array_to_view(pdel_in.data(),      ncol, nlev,  pdel_2d);
  array_to_view(inv_exner_in.data(), ncol, nlev,  inv_exner_2d);
  array_to_view(thv_in.data(),       ncol, nlev,  thv_2d);

  SHOC::SHOCInput shoc_input{host_dx_1d, host_dy_1d, zt_grid_2d, zi_grid_2d,
                            pres_2d, presi_2d, pdel_2d, thv_2d,
                            w_field_2d, wthl_sfc_1d, wqw_sfc_1d, uw_sfc_1d,
                            vw_sfc_1d, wtracer_sfc_2d, inv_exner_2d, phis_1d};

  reshape(host_dse,     host_dse_in.data(),     ncol, nlev);
  reshape(tke,          tke_in.data(),          ncol, nlev);
  reshape(thetal,       thetal_in.data(),       ncol, nlev);
  reshape(qw,           qw_in.data(),           ncol, nlev);
  reshape(shoc_ql,      shoc_ql_in.data(),      ncol, nlev);
  reshape(wthv_sec,     wthv_sec_in.data(),     ncol, nlev);
  reshape(tk,           tk_in.data(),           ncol, nlev);
  reshape(tkh,          tkh_in.data(),          ncol, nlev);
  reshape(shoc_cldfrac, shoc_cldfrac_in.data(), ncol, nlev);

  view_2d host_dse_2d("host_dse", ncol, npack);                 // dry static energy [J/kg] : dse = Cp*T + g*z + phis
  view_2d tke_2d("tke", ncol, npack);                           // turbulent kinetic energy [m2/s2]
  view_2d thetal_2d("thetal", ncol, npack);                     // liquid water potential temperature [K]
  view_2d qw_2d("qw", ncol, npack);                             // total water mixing ratio [kg/kg]
  view_2d shoc_ql_2d("shoc_ql", ncol, npack);                   // Vector-valued wind (u,v) [m/s]
  view_2d wthv_sec_2d("wthv_sec", ncol, npack);                 // buoyancy flux [K m/s]
  view_2d tk_2d("tk", ncol, npack);                             // tracers [varies]
  view_2d tkh_2d("tkh", ncol, npack);                           // eddy coefficient for momentum [m2/s]
  view_2d shoc_cldfrac_2d("shoc_cldfrac", ncol, npack);         // eddy heat conductivity
  view_3d shoc_hwind_3d("shoc_hwind",ncol,2,npack);             // Cloud fraction [-]
  view_3d qtracers_3d("qtracers",ncol,num_qtracers,npack);       // cloud liquid mixing ratio [kg/kg]

  array_to_view(shoc_dse_in.data(),     ncol, nlev, host_dse_2d);
  array_to_view(tke_in.data(),          ncol, nlev, tke_2d);
  array_to_view(qw_in.data(),           ncol, nlev, qw_2d);
  array_to_view(shoc_ql_in.data(),      ncol, nlev, shoc_ql_2d);
  array_to_view(tk_in.data(),           ncol, nlev, tk_2d);
  array_to_view(tkh_in.data(),          ncol, nlev, tkh_2d);
  array_to_view(wthv_sec_in.data(),     ncol, nlev, wthv_sec_2d);
  array_to_view(shoc_cldfrac_in.data(), ncol, nlev, shoc_cldfrac_2d);
  array_to_view(shoc_hwind_in.data(),   ncol, 2, nlev, shoc_hwind_3d);
  array_to_view(qtracers_in.data(),     ncol, num_qtracers, nlev, qtracers_3d);

  SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, qw_2d,
                                         shoc_hwind_3d, wthv_sec_2d, qtracers_3d,
                                         tk_2d, shoc_cldfrac_2d, shoc_ql_2d};

  view_1d pblh_1d("pblh",ncol);
  view_2d shoc_ql2_2d("shoc_ql2",ncol,npack);
  SHOC::SHOCOutput shoc_output{pblh_1d, shoc_ql2_2d};

  view_2d shoc_mix_2d("shoc_mix", ncol, npack),     // Turbulent length scale [m]
          w_sec_2d("w_sec", ncol, npack),           // vertical velocity variance [m2/s2]
          thl_sec_2d("thl_sec", ncol, nipack),      // temperature variance [K^2]
          qw_sec_2d("qw_sec", ncol, nipack),        // moisture variance [kg2/kg2]
          qwthl_sec_2d("qwthl_sec", ncol, nipack),  // temp moisture covariance [K kg/kg]
          wthl_sec_2d("wthl_sec", ncol, nipack),    // vertical heat flux [K m/s]
          wqw_sec_2d("wqw_sec", ncol, nipack),      // vertical moisture flux [K m/s]
          wtke_sec_2d("wtke_sec", ncol, nipack),    // vertical tke flux [m3/s3]
          uw_sec_2d("uw_sec", ncol, nipack),        // vertical zonal momentum flux [m2/s2]
          vw_sec_2d("vw_sec", ncol, nipack),        // vertical meridional momentum flux [m2/s2]
          w3_2d("w3", ncol, nipack),                // third moment vertical velocity [m3/s3]
          wqls_sec_2d("wqls_sec", ncol, npack),     // liquid water flux [kg/kg m/s]
          brunt_2d("brunt", ncol, npack),           // brunt vaisala frequency [s-1]
          isotropy_2d("isotropy", ncol, npack);     // return to isotropic timescale [s]

  SHOC::SHOCHistoryOutput shoc_history_output{shoc_mix_2d, w_sec_2d, thl_sec_2d, qw_sec_2d,
                                              qwthl_sec_2d, wthl_sec_2d, wqw_sec_2d, wtke_sec_2d,
                                              uw_sec_2d, vw_sec_2d, w3_2d, wqls_sec_2d, brunt_2d, isotropy_2d};

  const int nwind = ekat::npack<Spack>(2)*Spack::n;
  const int ntrac = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
  const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
  ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 128+(nwind+ntrac), policy);

  const auto elapsed_microsec = SHOC::shoc_main(ncol, nlev, nlevi, nlev, 1, num_qtracers, dtime, workspace_mgr,
                                                shoc_input, shoc_input_output, shoc_output, shoc_history_output);

  // get SHOC output back to CRM 
  view_to_array(shoc_input_output.tk,   ncol, nlev, tk_in);
  //view_to_array(shoc_input_output.tkh,  ncol, nlev, tkh_in);
  view_to_array(shoc_input_output.wthv_sec, ncol, nlev, wthv_sec_in);
  // view_to_array(shoc_input_output.tke, ncol, nlev, tke);
  view_to_array(shoc_input_output.shoc_cldfrac, ncol, nlev, shoc_cldfrac_in);
  view_to_array(shoc_input_output.horiz_wind, ncol, 2, nlev, shoc_hwind_in);
  view_to_array(shoc_input_output.host_dse, ncol, nlev, shoc_dse_in);
  view_to_array(shoc_input_output.qtracers, ncol, num_qtracers, nlev, qtracers_in);

  reshape(tk_in.data(),           tk,           ncol, nlev);
  reshape(wthv_sec_in.data(),     wthv_sec,     ncol, nlev);
  reshape(shoc_cldfrac_in.data(), shoc_cldfrac, ncol, nlev);
  reshape(shoc_dse_in.data(),     host_dse,     ncol, nlev);

  // update tracers and TKE
  parallel_for( SimpleBounds<3>(ncol,num_qtracers,nlev) , YAKL_LAMBDA (int i, int q, int k) {
    const int offset = (i*nlev+k)*num_qtracers+q;
    qtracers[offset] = qtracers_in(i,q,k);
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlev}), KOKKOS_LAMBDA(int i, int k) {
     const int offset = i*nlev+k;
     u_wind[offset] = shoc_hwind_in(i,0,k);
     v_wind[offset] = shoc_hwind_in(i,1,k);
  });
}
