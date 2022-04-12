
#pragma once

#include "pam_coupler.h"
#include "idealized_profiles.h"

inline void supercell_init( pam::PamCoupler &coupler ) {
  using yakl::c::parallel_for;
  using yakl::c::SimpleBounds;
  using profiles::init_supercell_temperature ;
  using profiles::init_supercell_pressure_dry;
  using profiles::init_supercell_sat_mix_dry ;
  using profiles::init_supercell_relhum      ;

  int constexpr ord = 5;
  SArray<real,1,ord> gll_pts, gll_wts;

  gll_pts(0)=-0.50000000000000000000000000000000000000;
  gll_pts(1)=-0.32732683535398857189914622812342917778;
  gll_pts(2)=0.00000000000000000000000000000000000000;
  gll_pts(3)=0.32732683535398857189914622812342917778;
  gll_pts(4)=0.50000000000000000000000000000000000000;

  gll_wts(0)=0.050000000000000000000000000000000000000;
  gll_wts(1)=0.27222222222222222222222222222222222222;
  gll_wts(2)=0.35555555555555555555555555555555555556;
  gll_wts(3)=0.27222222222222222222222222222222222222;
  gll_wts(4)=0.050000000000000000000000000000000000000;

  auto Rd    = coupler.get_R_d ();
  auto cp    = coupler.get_cp_d();
  auto p0    = coupler.get_p0  ();
  auto Rv    = coupler.get_R_v ();
  auto grav  = coupler.get_grav();

  auto &dm = coupler.get_data_manager_readwrite();

  auto vert_interface = dm.get<real,2>("vertical_interface_height");

  // This uses a piecewise linear profile for Temperature
  real constexpr z_0    = 0;
  real constexpr z_trop = 12000;
  real constexpr T_0    = 300;
  real constexpr T_trop = 213;
  real constexpr T_top  = 213;
  real constexpr p_0    = 100000;

  real4d quad_temp     ("quad_temp"     ,nz,ord-1,ord,nens);
  real3d hyPressureGLL ("hyPressureGLL" ,nz,ord,nens);

  real3d hyDensVapGLL  ("hyDensVapGLL"  ,nz,ord,nens);
  real2d hyDensVapCells("hyDensVapCells",nz,nens);

  // We'll use T and qv plus hydrostasis and equation of state to integrate total pressure
  // Combine  dp/dz=-rho_d*(1+qv)*g  with  p=rho_d*(Rd + qv*Rv)*T  via rho_d to create quadrature term for
  //     integrating ln(p)
  parallel_for( SimpleBounds<4>(nz,ord-1,ord,nens) ,
                YAKL_LAMBDA (int k, int kk, int kkk, int iens) {
    real dz = vert_interface(k+1,iens) - vert_interface(k,iens);
    // Middle of this cell
    real cellmid   = vert_interface(k,iens) + 0.5_fp*dz;
    // Bottom, top, and middle of the space between these two ord GLL points
    real ord_b    = cellmid + gll_pts(kk  )*dz;
    real ord_t    = cellmid + gll_pts(kk+1)*dz;
    real ord_m    = 0.5_fp * (ord_b + ord_t);
    // Compute grid spacing between these ord GLL points
    real ord_dz   = dz * ( gll_pts(kk+1) - gll_pts(kk) );
    // Compute the locate of this ord GLL point within the ord GLL points
    real zloc      = ord_m + ord_dz * gll_pts(kkk);
    // Compute full density at this location
    real ztop = vert_interface(nz,iens);
    real temp      = init_supercell_temperature (zloc, z_0, z_trop, ztop, T_0, T_trop, T_top);
    real press_dry = init_supercell_pressure_dry(zloc, z_0, z_trop, ztop, T_0, T_trop, T_top, p_0, Rd, grav);
    real qvs       = init_supercell_sat_mix_dry(press_dry, temp);
    real relhum    = init_supercell_relhum(zloc, z_0, z_trop);
    if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
    real qv        = min( 0.014_fp , qvs*relhum );
    quad_temp(k,kk,kkk,iens) = -(1+qv)*grav/(Rd+qv*Rv)/temp;
  });

  // Integrate to get total pressure at GLL points within each cell domain
  parallel_for( nens , YAKL_LAMBDA (int iens) {
    hyPressureGLL(0,0,iens) = p_0;
    for (int k=0; k < nz; k++) {
      for (int kk=0; kk < ord-1; kk++) {
        real tot = 0;
        for (int kkk=0; kkk < ord; kkk++) {
          tot += quad_temp(k,kk,kkk,iens) * gll_wts(kkk);
        }
        tot *= dz(k,iens) * ( gll_pts(kk+1) - gll_pts(kk) );
        hyPressureGLL(k,kk+1,iens) = hyPressureGLL(k,kk,iens) * exp( tot );
        if (kk == ord-2 && k < nz-1) {
          hyPressureGLL(k+1,0,iens) = hyPressureGLL(k,ord-1,iens);
        }
      }
    }
  });

  real2d rho_d_col("rho_d_col",nz,nens);
  real2d uvel_col ("uvel_col" ,nz,nens);
  real2d vvel_col ("vvel_col" ,nz,nens);
  real2d wvel_col ("wvel_col" ,nz,nens);
  real2d temp_col ("temp_col" ,nz,nens);
  real2d rho_v_col("rho_v_col",nz,nens);

  // Compute cell averages of dry density, wind velocities, temperature, and water vapor mass
  parallel_for( SimpleBounds<2>(nz,nens) , YAKL_LAMBDA (int k, int iens) {
    rho_d_col = 0;
    uvel_col  = 0;
    vvel_col  = 0;
    wvel_col  = 0;
    temp_col  = 0;
    rho_v_col = 0;
    for (int kk=0; kk < ord; kk++) {
      real dz        = vert_interface(k+1,iens) - vert_interface(k,iens);
      real zmid      = 0.5_fp * (vert_interface(k,iens) + vert_interface(k+1,iens));
      real zloc      = zmid + gll_pts(kk)*dz;
      real ztop      = vert_interface(nz,iens);
      real temp      = init_supercell_temperature (zloc, z_0, z_trop, ztop, T_0, T_trop, T_top);
      real press_dry = init_supercell_pressure_dry(zloc, z_0, z_trop, ztop, T_0, T_trop, T_top, p_0, Rd, grav);
      real qvs       = init_supercell_sat_mix_dry(press_dry, temp);
      real relhum    = init_supercell_relhum(zloc, z_0, z_trop);
      if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
      real qv        = min( 0.014_fp , qvs*relhum );
      real p         = hyPressureGLL(k,kk,iens);
      real rho_d     = p / (Rd + qv*Rv) / temp;
      real rho_v     = qv * rho_d;

      real uvel;
      real constexpr zs = 5000;
      real constexpr us = 30;
      real constexpr uc = 15;
      if (zloc < zs) { uvel = us * (zloc / zs) - uc; }
      else           { uvel = us - uc;               }
      real vvel       = 0;
      real wvel       = 0;

      rho_d_col(k,iens) += rho_d * gll_wts(kk);
      uvel_col (k,iens) += uvel  * gll_wts(kk);
      vvel_col (k,iens) += vvel  * gll_wts(kk);
      wvel_col (k,iens) += wvel  * gll_wts(kk);
      temp_col (k,iens) += temp  * gll_wts(kk);
      rho_v_col(k,iens) += rho_v * gll_wts(kk);
    }
  });

  auto rho_d = dm.get<real,4>("density_dry");
  auto uvel  = dm.get<real,4>("uvel"       );
  auto vvel  = dm.get<real,4>("vvel"       );
  auto wvel  = dm.get<real,4>("wvel"       );
  auto temp  = dm.get<real,4>("temp"       );
  auto rho_v = dm.get<real,4>("water_vapor");

  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
    rho_d(k,j,i,iens) = rho_d_col(k,iens);
    uvel (k,j,i,iens) = uvel_col (k,iens);
    vvel (k,j,i,iens) = vvel_col (k,iens);
    wvel (k,j,i,iens) = wvel_col (k,iens);
    temp (k,j,i,iens) = temp_col (k,iens);
    rho_v(k,j,i,iens) = rho_v_col(k,iens);
  });

}


