
#pragma once

#include "pam_utils.h"

#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "shoc_f90.hpp"

using namespace scream;
using namespace scream::shoc;

using array1D = IRType<Real, 1>::type;
using array2D = IRType<Real, 2>::type;
using array3D = IRType<Real, 3>::type;

void shoc_main_cxx(int &shcol, int &nlev, int &nlevi, double &dtime, int &nadv, array1D& host_dx, array1D& host_dy,
                   array2D& thv, array2D& zt_grid, array2D& zi_grid, array2D& pres, array2D& presi, array2D& pdel,
                   array1D& wthl_sfc, array1D& wqw_sfc, array1D& uw_sfc, array1D& vw_sfc, array2D& wtracer_sfc,
                   int &num_qtracers, array2D& w_field, array2D& inv_exner, array1D& phis, array2D& host_dse, array2D& tke,
                   array2D& thetal, array2D& qw, array2D& u_wind, array2D& v_wind, array3D& qtracers,array2D& wthv_sec,
                   array2D& tkh, array2D& tk, array2D& shoc_ql, array2D& shoc_cldfrac, array1D& pblh, array2D& shoc_mix,
                   array2D& isotropy, array2D& w_sec, array2D& thl_sec, array2D& qw_sec, array2D& qwthl_sec,
                   array2D& wthl_sec, array2D& wqw_sec, array2D& wtke_sec, array2D& uw_sec, array2D& vw_sec, array2D& w3,
                   array2D& wqls_sec, array2D& brunt, array2D& shoc_ql2 )

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

  // -------------------------------------------------
  // Set surface geopotential and fluxes
  // -------------------------------------------------
  view_1d host_dx_1d("host_dx", ncol),
          host_dy_1d("host_dy", ncol),
          wthl_sfc_1d("wthl_sfc", ncol),  // Surface sensible heat flux [K m/s]
          wqw_sfc_1d("wqw_sfc", ncol),    // Surface latent heat flux [kg/kg m/s]
          uw_sfc_1d("uw_sfc", ncol),      // Surface momentum flux (u-direction) [m2/s2]
          vw_sfc_1d("vw_sfc", ncol),      // Surface momentum flux (v-direction) [m2/s2]
          phis_1d("phis_1d", ncol);       // Host model surface geopotential height

  Kokkos::parallel_for("host grid", ncol, KOKKOS_LAMBDA (const int& i) {
    host_dx_1d(i)  = host_dx.data()[i];
    host_dy_1d(i)  = host_dy.data()[i];
    wthl_sfc_1d(i) = wthl_sfc.data()[i];
    wqw_sfc_1d(i)  = wqw_sfc.data()[i];
    uw_sfc_1d(i)   = uw_sfc.data()[i];
    vw_sfc_1d(i)   = vw_sfc.data()[i];
    phis_1d(i)     = phis.data()[i];
  });

  auto zt_g_in        = reshape(zt_grid);
  auto zi_g_in        = reshape(zi_grid);
  auto pres_in        = reshape(pres);
  auto presi_in       = reshape(presi);
  auto pdel_in        = reshape(pdel);
  auto inv_exner_in   = reshape(inv_exner);
  auto thv_in         = reshape(thv);
  auto w_field_in     = reshape(w_field);
  auto wtracer_sfc_in = reshape(wtracer_sfc);

  view_2d zt_grid_2d("zt_grid", ncol, npack),                 // heights, for thermo grid [m]
          zi_grid_2d("zi_grid", ncol, nipack),                // heights, for interface grid [m]
          pres_2d("pres", ncol, npack),                       // pressure levels on thermo grid [Pa]
          presi_2d("presi", ncol, nipack),                    // pressure levels on interface grid [Pa]
          pdel_2d("pdel", ncol, npack),                       // Differences in pressure levels [Pa]
          thv_2d("thv", ncol, npack),                         // virtual potential temperature [K]
          w_field_2d("w_field", ncol, npack),                 // large scale vertical velocity [m/s]
          wtracer_sfc_2d("wtracer", ncol, num_qtracers),             // Surface flux for tracers [varies]
          inv_exner_2d("inv_exner", ncol, npack);

  array_to_view(zt_g_in.data(),        ncol, nlev,  zt_grid_2d);
  array_to_view(zi_g_in.data(),        ncol, nlevi, zi_grid_2d);
  array_to_view(pres_in.data(),        ncol, nlev,  pres_2d);
  array_to_view(presi_in.data(),       ncol, nlevi, presi_2d);
  array_to_view(pdel_in.data(),        ncol, nlev,  pdel_2d);
  array_to_view(inv_exner_in.data(),   ncol, nlev,  inv_exner_2d);
  array_to_view(thv_in.data(),         ncol, nlev,  thv_2d);
  array_to_view(wtracer_sfc_in.data(), ncol, num_qtracers,  wtracer_sfc_2d);
  array_to_view(w_field_in.data(),     ncol, nlev,  w_field_2d);

  SHOC::SHOCInput shoc_input{host_dx_1d, host_dy_1d, zt_grid_2d, zi_grid_2d,
                            pres_2d, presi_2d, pdel_2d, thv_2d,
                            w_field_2d, wthl_sfc_1d, wqw_sfc_1d, uw_sfc_1d,
                            vw_sfc_1d, wtracer_sfc_2d, inv_exner_2d, phis_1d};

  auto host_dse_in     = reshape(host_dse);
  auto tke_in          = reshape(tke);
  auto thetal_in       = reshape(thetal);
  auto qw_in           = reshape(qw);
  auto shoc_ql_in      = reshape(shoc_ql);
  auto wthv_sec_in     = reshape(wthv_sec);
  auto tk_in           = reshape(tk);
  auto shoc_cldfrac_in = reshape(shoc_cldfrac);
  auto tkh_in          = reshape(tkh);

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

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int k = ilev*Spack::n + s;
    int offset = icol*nlev+k;
    if (k < nlev) {
     shoc_hwind_3d(icol,0,ilev)[s] = u_wind.data()[offset];
     shoc_hwind_3d(icol,1,ilev)[s] = v_wind.data()[offset];
    }
  });

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0}, {ncol, num_qtracers, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int q, int ilev, int s) {
    int k = ilev*Spack::n + s;
    const int offset = (icol*nlev+k)*num_qtracers+q;
    if (k < nlev) {
      qtracers_3d(icol,q,ilev)[s] = qtracers.data()[offset];
    }
  });

  array_to_view(host_dse_in.data(),     ncol, nlev, host_dse_2d);
  array_to_view(tke_in.data(),          ncol, nlev, tke_2d);
  array_to_view(thetal_in.data(),       ncol, nlev, thetal_2d);
  array_to_view(qw_in.data(),           ncol, nlev, qw_2d);
  array_to_view(shoc_ql_in.data(),      ncol, nlev, shoc_ql_2d);
  array_to_view(tk_in.data(),           ncol, nlev, tk_2d);
  array_to_view(tkh_in.data(),          ncol, nlev, tkh_2d);
  array_to_view(wthv_sec_in.data(),     ncol, nlev, wthv_sec_2d);
  array_to_view(shoc_cldfrac_in.data(), ncol, nlev, shoc_cldfrac_2d);
  array_to_view(shoc_ql_in.data(),      ncol, nlev, shoc_ql_2d);

  SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, qw_2d,
                                         shoc_hwind_3d, wthv_sec_2d, qtracers_3d,
                                         tk_2d, shoc_cldfrac_2d, shoc_ql_2d};

  view_1d pblh_1d("pblh",ncol);
  view_2d shoc_ql2_2d("shoc_ql2",ncol, npack);
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
  ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 13+(nwind+ntrac), policy);

  const auto elapsed_microsec = SHOC::shoc_main(shcol, nlev, nlevi, nlev, nadv, num_qtracers, dtime, workspace_mgr,
                                                shoc_input, shoc_input_output, shoc_output, shoc_history_output);

  // update tk, tke, and other outputs
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int k = ilev*Spack::n + s;
    int offset = icol*nlev+k;
    if (k < nlev) {
      tk.data()[offset]           = shoc_input_output.tk(icol,ilev)[s];
      tke.data()[offset]          = shoc_input_output.tke(icol,ilev)[s];
      wthv_sec.data()[offset]     = shoc_input_output.wthv_sec(icol,ilev)[s];
      host_dse.data()[offset]     = shoc_input_output.host_dse(icol,ilev)[s];
      shoc_cldfrac.data()[offset] = shoc_input_output.shoc_cldfrac(icol,ilev)[s];
      u_wind.data()[offset]       = shoc_input_output.horiz_wind(icol,0,ilev)[s];
      v_wind.data()[offset]       = shoc_input_output.horiz_wind(icol,1,ilev)[s];
    }
  });

  // update tracers
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0}, {ncol, num_qtracers, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int q, int ilev, int s) {
    int k = ilev*Spack::n + s;
    const int offset = (icol*nlev+k)*num_qtracers+q;
    if (k < nlev) {
      qtracers.data()[offset] = shoc_input_output.qtracers(icol,q,ilev)[s];
    }
  });
}
