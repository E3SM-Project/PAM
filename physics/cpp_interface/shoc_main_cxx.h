
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
          phis_1d("phis_1d", ncol),       // Host model surface geopotential height
          pblh_1d("pblh", ncol);


  Kokkos::parallel_for("host grid", ncol, KOKKOS_LAMBDA (const int& i) {
    host_dx_1d(i)  = host_dx.data()[i];
    host_dy_1d(i)  = host_dy.data()[i];
    wthl_sfc_1d(i) = wthl_sfc.data()[i];
    wqw_sfc_1d(i)  = wqw_sfc.data()[i];
    uw_sfc_1d(i)   = uw_sfc.data()[i];
    vw_sfc_1d(i)   = vw_sfc.data()[i];
    phis_1d(i)     = phis.data()[i];
    pblh_1d(i)     = pblh.data()[i];
  });

  view_2d zt_grid_2d("zt_grid", ncol, npack),                 // heights, for thermo grid [m]
          zi_grid_2d("zi_grid", ncol, nipack),                // heights, for interface grid [m]
          pres_2d("pres", ncol, npack),                       // pressure levels on thermo grid [Pa]
          presi_2d("presi", ncol, nipack),                    // pressure levels on interface grid [Pa]
          pdel_2d("pdel", ncol, npack),                       // Differences in pressure levels [Pa]
          thv_2d("thv", ncol, npack),                         // virtual potential temperature [K]
          w_field_2d("w_field", ncol, npack),                 // large scale vertical velocity [m/s]
          wtracer_sfc_2d("wtracer", ncol, num_qtracers),             // Surface flux for tracers [varies]
          inv_exner_2d("inv_exner", ncol, npack);

  array_to_view(zt_grid.data(),     DataFormat::NCWH, ncol, nlev,  zt_grid_2d);
  array_to_view(zi_grid.data(),     DataFormat::NCWH, ncol, nlevi, zi_grid_2d);
  array_to_view(pres.data(),        DataFormat::NCWH, ncol, nlev,  pres_2d);
  array_to_view(presi.data(),       DataFormat::NCWH, ncol, nlevi, presi_2d);
  array_to_view(pdel.data(),        DataFormat::NCWH, ncol, nlev,  pdel_2d);
  array_to_view(inv_exner.data(),   DataFormat::NCWH, ncol, nlev,  inv_exner_2d);
  array_to_view(thv.data(),         DataFormat::NCWH, ncol, nlev,  thv_2d);
  array_to_view(wtracer_sfc.data(), DataFormat::NCWH, ncol, num_qtracers, wtracer_sfc_2d);
  array_to_view(w_field.data(),     DataFormat::NCWH, ncol, nlev,  w_field_2d);

//  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
//    int off = (ilev*Spack::n+s)*ncol+icol;
//    printf("%d, %d, %13.6f, %13.6f\n",icol,ilev,zt_grid_2d(icol,ilev)[s],zt_grid.data()[off ]);
//  });

  SHOC::SHOCInput shoc_input{host_dx_1d, host_dy_1d, zt_grid_2d, zi_grid_2d,
                            pres_2d, presi_2d, pdel_2d, thv_2d,
                            w_field_2d, wthl_sfc_1d, wqw_sfc_1d, uw_sfc_1d,
                            vw_sfc_1d, wtracer_sfc_2d, inv_exner_2d, phis_1d};

  view_2d host_dse_2d("host_dse", ncol, npack);                 // dry static energy [J/kg] : dse = Cp*T + g*z + phis
  view_2d tke_2d("tke", ncol, npack);                           // turbulent kinetic energy [m2/s2]
  view_2d thetal_2d("thetal", ncol, npack);                     // liquid water potential temperature [K]
  view_2d qw_2d("qw", ncol, npack);                             // total water mixing ratio [kg/kg]
  view_2d shoc_ql_2d("shoc_ql", ncol, npack);                   // Vector-valued wind (u,v) [m/s]
  view_2d wthv_sec_2d("wthv_sec", ncol, npack);                 // buoyancy flux [K m/s]
  view_2d tk_2d("tk", ncol, npack);                             // tracers [varies]
  view_2d tkh_2d("tkh", ncol, npack);                           // eddy coefficient for momentum [m2/s]
  view_2d shoc_cldfrac_2d("shoc_cldfrac", ncol, npack);         // eddy heat conductivity
  view_2d u_wind_2d("u_wind", ncol, npack);
  view_2d v_wind_2d("v_wind", ncol, npack);
  view_3d shoc_hwind_3d("shoc_hwind",ncol,2,npack);             // Cloud fraction [-]
  view_3d qtracers_3d("qtracers",ncol,num_qtracers,npack);      // cloud liquid mixing ratio [kg/kg]

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
    int k = ilev*Spack::n + s;
    int offset = k*ncol + icol;
    if (k < nlev) {
     shoc_hwind_3d(icol,0,ilev)[s] = u_wind.data()[offset];
     shoc_hwind_3d(icol,1,ilev)[s] = v_wind.data()[offset];
    }
  });

#if 0
  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0}, {ncol, num_qtracers, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int q, int ilev, int s) {
    int k = ilev*Spack::n + s;
    const int offset = (q*nlev+k)*ncol+icol;
    if (k < nlev) {
      qtracers_3d(icol,q,ilev)[s] = qtracers.data()[offset];
    }
  });
#endif

  array_to_view(host_dse.data(),     DataFormat::NCWH, ncol, nlev, host_dse_2d);
  array_to_view(tke.data(),          DataFormat::NCWH, ncol, nlev, tke_2d);
  array_to_view(thetal.data(),       DataFormat::NCWH, ncol, nlev, thetal_2d);
  array_to_view(qw.data(),           DataFormat::NCWH, ncol, nlev, qw_2d);
  array_to_view(shoc_ql.data(),      DataFormat::NCWH, ncol, nlev, shoc_ql_2d);
  array_to_view(tk.data(),           DataFormat::NCWH, ncol, nlev, tk_2d);
  array_to_view(tkh.data(),          DataFormat::NCWH, ncol, nlev, tkh_2d);
  array_to_view(wthv_sec.data(),     DataFormat::NCWH, ncol, nlev, wthv_sec_2d);
  array_to_view(shoc_cldfrac.data(), DataFormat::NCWH, ncol, nlev, shoc_cldfrac_2d);
  array_to_view(shoc_ql.data(),      DataFormat::NCWH, ncol, nlev, shoc_ql_2d);
  array_to_view(u_wind.data(),       DataFormat::NCWH, ncol, nlev, u_wind_2d);
  array_to_view(v_wind.data(),       DataFormat::NCWH, ncol, nlev, v_wind_2d);
  array_to_view(qtracers.data(),     DataFormat::NWCH, ncol, num_qtracers, nlev, qtracers_3d);

  Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {ncol, npack, Spack::n}), KOKKOS_LAMBDA(int icol, int ilev, int s) {
     shoc_hwind_3d(icol,0,ilev)[s] = u_wind_2d(icol, ilev)[s];
     shoc_hwind_3d(icol,1,ilev)[s] = v_wind_2d(icol, ilev)[s];
  });

  SHOC::SHOCInputOutput shoc_input_output{host_dse_2d, tke_2d, thetal_2d, qw_2d,
                                         shoc_hwind_3d, wthv_sec_2d, qtracers_3d,
                                         tk_2d, shoc_cldfrac_2d, shoc_ql_2d};

  view_2d shoc_ql2_2d("shoc_ql2",ncol, npack);
  array_to_view(shoc_ql2.data(), DataFormat::NCWH, ncol, nlev, shoc_ql2_2d);

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

  array_to_view(shoc_mix.data(),  DataFormat::NCWH, ncol, nlev, shoc_mix_2d);
  array_to_view(w_sec.data(),     DataFormat::NCWH, ncol, nlev, w_sec_2d);
  array_to_view(thl_sec.data(),   DataFormat::NCWH, ncol, nlevi, thl_sec_2d);
  array_to_view(qw_sec.data(),    DataFormat::NCWH, ncol, nlevi, qw_sec_2d);
  array_to_view(qwthl_sec.data(), DataFormat::NCWH, ncol, nlevi, qwthl_sec_2d);
  array_to_view(wthl_sec.data(),  DataFormat::NCWH, ncol, nlevi, wthl_sec_2d);
  array_to_view(wqw_sec.data(),   DataFormat::NCWH, ncol, nlevi, wqw_sec_2d);
  array_to_view(wtke_sec.data(),  DataFormat::NCWH, ncol, nlevi, wtke_sec_2d);
  array_to_view(uw_sec.data(),    DataFormat::NCWH, ncol, nlevi, uw_sec_2d);
  array_to_view(vw_sec.data(),    DataFormat::NCWH, ncol, nlevi, vw_sec_2d);
  array_to_view(w3.data(),        DataFormat::NCWH, ncol, nlevi, w3_2d);
  array_to_view(wqls_sec.data(),  DataFormat::NCWH, ncol, nlev, wqls_sec_2d);
  array_to_view(brunt.data(),     DataFormat::NCWH, ncol, nlev, brunt_2d);
  array_to_view(isotropy.data(),  DataFormat::NCWH, ncol, nlev, isotropy_2d);

  SHOC::SHOCHistoryOutput shoc_history_output{shoc_mix_2d, w_sec_2d, thl_sec_2d, qw_sec_2d,
                                              qwthl_sec_2d, wthl_sec_2d, wqw_sec_2d, wtke_sec_2d,
                                              uw_sec_2d, vw_sec_2d, w3_2d, wqls_sec_2d, brunt_2d, isotropy_2d};

  const int nwind = ekat::npack<Spack>(2)*Spack::n;
  const int ntrac = ekat::npack<Spack>(num_qtracers+3)*Spack::n;
  const auto policy = ekat::ExeSpaceUtils<SHOC::KT::ExeSpace>::get_default_team_policy(ncol, npack);
  ekat::WorkspaceManager<Spack, SHOC::KT::Device> workspace_mgr(nipack, 13+(nwind+ntrac), policy);

  const auto elapsed_microsec = SHOC::shoc_main(shcol, nlev, nlevi, nlev, nadv, num_qtracers, dtime, workspace_mgr,
                                                shoc_input, shoc_input_output, shoc_output, shoc_history_output);

  // get SHOC output back to CRM 
  view_to_array(shoc_input_output.tk,           ncol, nlev, DataFormat::NCWH, tk);
  view_to_array(shoc_input_output.wthv_sec,     ncol, nlev, DataFormat::NCWH, wthv_sec);
  view_to_array(shoc_input_output.tke,          ncol, nlev, DataFormat::NCWH, tke);
  view_to_array(shoc_input_output.shoc_cldfrac, ncol, nlev, DataFormat::NCWH, shoc_cldfrac);
//  view_to_array(shoc_input_output.horiz_wind, ncol, 2, nlev, DataFormat::NCWH, shoc_hwind);
  view_to_array(shoc_input_output.host_dse,     ncol, nlev, DataFormat::NCWH, host_dse);
  view_to_array(shoc_input_output.qtracers,     ncol, num_qtracers, nlev, DataFormat::NWCH, qtracers);
  view_to_array(shoc_input_output.qw,           ncol, nlev, DataFormat::NCWH, qw);
  view_to_array(shoc_input_output.shoc_ql,      ncol, nlev, DataFormat::NCWH, shoc_ql);
  view_to_array(shoc_input_output.thetal,       ncol, nlev, DataFormat::NCWH, thetal);
  view_to_array(shoc_output.shoc_ql2,           ncol, nlev, DataFormat::NCWH, shoc_ql2);
}
