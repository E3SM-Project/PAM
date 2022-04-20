
#pragma once

#include "pam_const.h"
#include "idealized_profiles.h"

inline void supercell_init( realConst1d vert_interface , real1d &rho_d_col , real1d &uvel_col , real1d &vvel_col ,
                            real1d &wvel_col , real1d &temp_col , real1d &rho_v_col , real Rd ,
                            real Rv , real grav) {
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

  // This uses a piecewise linear profile for Temperature
  real constexpr z_0    = 0;
  real constexpr z_trop = 12000;
  real constexpr T_0    = 300;
  real constexpr T_trop = 213;
  real constexpr T_top  = 213;
  real constexpr p_0    = 100000;

  int nz = yakl::intrinsics::size(vert_interface)-1;

  real3d quad_temp("quad_temp",nz,ord-1,ord);

  // We'll use T and qv plus hydrostasis and equation of state to integrate total pressure
  // Combine  dp/dz=-rho_d*(1+qv)*g  with  p=rho_d*(Rd + qv*Rv)*T  via rho_d to create quadrature term for
  //     integrating ln(p)
  parallel_for( SimpleBounds<3>(nz,ord-1,ord) ,
                YAKL_LAMBDA (int k, int kk, int kkk) {
    real dz = vert_interface(k+1) - vert_interface(k);
    // Middle of this cell
    real cellmid   = vert_interface(k) + 0.5_fp*dz;
    // Bottom, top, and middle of the space between these two ord GLL points
    real ord_b    = cellmid + gll_pts(kk  )*dz;
    real ord_t    = cellmid + gll_pts(kk+1)*dz;
    real ord_m    = 0.5_fp * (ord_b + ord_t);
    // Compute grid spacing between these ord GLL points
    real ord_dz   = dz * ( gll_pts(kk+1) - gll_pts(kk) );
    // Compute the locate of this ord GLL point within the ord GLL points
    real zloc      = ord_m + ord_dz * gll_pts(kkk);
    // Compute full density at this location
    real ztop = vert_interface(nz);
    real temp      = init_supercell_temperature (zloc, z_0, z_trop, ztop, T_0, T_trop, T_top);
    real press_dry = init_supercell_pressure_dry(zloc, z_0, z_trop, ztop, T_0, T_trop, T_top, p_0, Rd, grav);
    real qvs       = init_supercell_sat_mix_dry(press_dry, temp);
    real relhum    = init_supercell_relhum(zloc, z_0, z_trop);
    if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
    real qv        = std::min( 0.014_fp , qvs*relhum );
    quad_temp(k,kk,kkk) = -(1+qv)*grav/(Rd+qv*Rv)/temp;
  });

  real2d hyPressureGLL("hyPressureGLL",nz,ord);

  // Integrate to get total pressure at GLL points within each cell domain
  parallel_for( 1 , YAKL_LAMBDA (int dummy) {
    hyPressureGLL(0,0) = p_0;
    for (int k=0; k < nz; k++) {
      real dz = vert_interface(k+1) - vert_interface(k);
      for (int kk=0; kk < ord-1; kk++) {
        real tot = 0;
        for (int kkk=0; kkk < ord; kkk++) {
          tot += quad_temp(k,kk,kkk) * gll_wts(kkk);
        }
        tot *= dz * ( gll_pts(kk+1) - gll_pts(kk) );
        hyPressureGLL(k,kk+1) = hyPressureGLL(k,kk) * exp( tot );
        if (kk == ord-2 && k < nz-1) {
          hyPressureGLL(k+1,0) = hyPressureGLL(k,ord-1);
        }
      }
    }
  });

  quad_temp = real3d();

  // Compute cell averages of dry density, wind velocities, temperature, and water vapor mass
  parallel_for( SimpleBounds<1>(nz) , YAKL_LAMBDA (int k) {
    rho_d_col(k) = 0;
    uvel_col (k) = 0;
    vvel_col (k) = 0;
    wvel_col (k) = 0;
    temp_col (k) = 0;
    rho_v_col(k) = 0;
    for (int kk=0; kk < ord; kk++) {
      real dz        = vert_interface(k+1) - vert_interface(k);
      real zmid      = 0.5_fp * (vert_interface(k) + vert_interface(k+1));
      real zloc      = zmid + gll_pts(kk)*dz;
      real ztop      = vert_interface(nz);
      real temp      = init_supercell_temperature (zloc, z_0, z_trop, ztop, T_0, T_trop, T_top);
      real press_dry = init_supercell_pressure_dry(zloc, z_0, z_trop, ztop, T_0, T_trop, T_top, p_0, Rd, grav);
      real qvs       = init_supercell_sat_mix_dry(press_dry, temp);
      real relhum    = init_supercell_relhum(zloc, z_0, z_trop);
      if (relhum * qvs > 0.014_fp) relhum = 0.014_fp / qvs;
      real qv        = std::min( 0.014_fp , qvs*relhum );
      real p         = hyPressureGLL(k,kk);
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

      rho_d_col(k) += rho_d * gll_wts(kk);
      uvel_col (k) += uvel  * gll_wts(kk);
      vvel_col (k) += vvel  * gll_wts(kk);
      wvel_col (k) += wvel  * gll_wts(kk);
      temp_col (k) += temp  * gll_wts(kk);
      rho_v_col(k) += rho_v * gll_wts(kk);
    }
  });

}


