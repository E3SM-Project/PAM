
#pragma once

#include "pam_const.h"
#include "DataManager.h"
#include "pam_coupler.h"


inline void compute_gcm_forcing_tendencies( PamCoupler &coupler , real2d &rho_d_in , real2d &uvel_in , real2d &vvel_in ,
                                            real2d &wvel_in , real2d &temp_in , real2d &rho_v_in , real dt_gcm ) {
  using yakl::atomicAdd;
  auto &dm = coupler.dm;
  // Get current state from coupler
  auto rho_d = dm.get<real const,4>( "density_dry" );
  auto uvel  = dm.get<real const,4>( "uvel"        );
  auto vvel  = dm.get<real const,4>( "vvel"        );
  auto wvel  = dm.get<real const,4>( "wvel"        );
  auto temp  = dm.get<real const,4>( "temp"        );
  auto rho_v = dm.get<real const,4>( "water_vapor" );

  int nz   = dm.get_dimension_size("z"   );
  int ny   = dm.get_dimension_size("y"   );
  int nx   = dm.get_dimension_size("x"   );
  int nens = dm.get_dimension_size("nens");

  real2d colavg_rho_d("colavg_rho_d",nz,nens);
  real2d colavg_uvel ("colavg_uvel" ,nz,nens);
  real2d colavg_vvel ("colavg_vvel" ,nz,nens);
  real2d colavg_wvel ("colavg_wvel" ,nz,nens);
  real2d colavg_temp ("colavg_temp" ,nz,nens);
  real2d colavg_rho_v("colavg_rho_v",nz,nens);

  parallel_for( "Initialize colsum to zero" , SimpleBounds<2>(nz,nens) , 
                YAKL_LAMBDA (int k, int iens) {
    colavg_rho_d(k,iens) = 0;
    colavg_uvel (k,iens) = 0;
    colavg_vvel (k,iens) = 0;
    colavg_wvel (k,iens) = 0;
    colavg_temp (k,iens) = 0;
    colavg_rho_v(k,iens) = 0;
  });

  real r_nx_ny  = 1._fp / (nx*ny);
  parallel_for( "Compute summed column of current state" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
    atomicAdd( colavg_rho_d(k,iens) , rho_d(k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_uvel (k,iens) , uvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_vvel (k,iens) , vvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_wvel (k,iens) , wvel (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_temp (k,iens) , temp (k,j,i,iens) * r_nx_ny );
    atomicAdd( colavg_rho_v(k,iens) , rho_v(k,j,i,iens) * r_nx_ny );
  });

  // We need the GCM forcing tendencies later, so store these in the coupler's data manager
  if (! dm.entry_exists("gcm_tend_rho_d")) {
    dm.register_and_allocate<real>( "gcm_tend_rho_d" , "GCM forcing tendency for dry density"         , {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_uvel"  , "GCM forcing tendency for u-velocity"          , {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_vvel"  , "GCM forcing tendency for v-velocity"          , {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_wvel"  , "GCM forcing tendency for w-velocity"          , {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_temp"  , "GCM forcing tendency for temperature"         , {nz,nens} , {"z","nens"} );
    dm.register_and_allocate<real>( "gcm_tend_rho_v" , "GCM forcing tendency for water vapor density" , {nz,nens} , {"z","nens"} );
  }

  // Retrieve the GCM forcing tendency arrays
  auto gcm_tend_rho_d = dm.get<real,2>("gcm_tend_rho_d");
  auto gcm_tend_uvel  = dm.get<real,2>("gcm_tend_uvel" );
  auto gcm_tend_vvel  = dm.get<real,2>("gcm_tend_vvel" );
  auto gcm_tend_wvel  = dm.get<real,2>("gcm_tend_wvel" );
  auto gcm_tend_temp  = dm.get<real,2>("gcm_tend_temp" );
  auto gcm_tend_rho_v = dm.get<real,2>("gcm_tend_rho_v");

  real r_dt_gcm = 1._fp / dt_gcm;
  parallel_for( "Compute GCM forcing tendencies" , SimpleBounds<2>(nz,nens) , 
                YAKL_LAMBDA (int k, int iens) {
    gcm_tend_rho_d(k,iens) = ( rho_d_in(k,iens) - colavg_rho_d(k,iens) ) * r_dt_gcm;
    gcm_tend_uvel (k,iens) = ( uvel_in (k,iens) - colavg_uvel (k,iens) ) * r_dt_gcm;
    gcm_tend_vvel (k,iens) = ( vvel_in (k,iens) - colavg_vvel (k,iens) ) * r_dt_gcm;
    gcm_tend_wvel (k,iens) = ( wvel_in (k,iens) - colavg_wvel (k,iens) ) * r_dt_gcm;
    gcm_tend_temp (k,iens) = ( temp_in (k,iens) - colavg_temp (k,iens) ) * r_dt_gcm;
    gcm_tend_rho_v(k,iens) = ( rho_v_in(k,iens) - colavg_rho_v(k,iens) ) * r_dt_gcm;
  });
}




inline void apply_gcm_forcing_tendencies( PamCoupler &coupler , real dt ) {
  using yakl::atomicAdd;
  auto &dm = coupler.dm;

  int nz   = dm.get_dimension_size("z"   );
  int ny   = dm.get_dimension_size("y"   );
  int nx   = dm.get_dimension_size("x"   );
  int nens = dm.get_dimension_size("nens");

  auto rho_d = dm.get<real,4>( "density_dry" );
  auto uvel  = dm.get<real,4>( "uvel"        );
  auto vvel  = dm.get<real,4>( "vvel"        );
  auto wvel  = dm.get<real,4>( "wvel"        );
  auto temp  = dm.get<real,4>( "temp"        );
  auto rho_v = dm.get<real,4>( "water_vapor" );

  auto gcm_tend_rho_d = dm.get<real,2>("gcm_tend_rho_d");
  auto gcm_tend_uvel  = dm.get<real,2>("gcm_tend_uvel" );
  auto gcm_tend_vvel  = dm.get<real,2>("gcm_tend_vvel" );
  auto gcm_tend_wvel  = dm.get<real,2>("gcm_tend_wvel" );
  auto gcm_tend_temp  = dm.get<real,2>("gcm_tend_temp" );
  auto gcm_tend_rho_v = dm.get<real,2>("gcm_tend_rho_v");

  real3d rho_v_neg_mass("rho_v_neg_mass",ny,nx,nens);
  real3d rho_v_pos_mass("rho_v_pos_mass",ny,nx,nens);

  parallel_for( Bounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int iens) {
    rho_v_neg_mass(j,i,iens) = 0;
    rho_v_pos_mass(j,i,iens) = 0;
  });

  parallel_for( "Apply GCM forcing" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_DEVICE_LAMBDA (int k, int j, int i, int iens) {
    // Apply forcing
    rho_d(k,j,i,iens) += gcm_tend_rho_d(k,iens) * dt;
    uvel (k,j,i,iens) += gcm_tend_uvel (k,iens) * dt;
    vvel (k,j,i,iens) += gcm_tend_vvel (k,iens) * dt;
    wvel (k,j,i,iens) += gcm_tend_wvel (k,iens) * dt;
    temp (k,j,i,iens) += gcm_tend_temp (k,iens) * dt;
    rho_v(k,j,i,iens) += gcm_tend_rho_v(k,iens) * dt;
    // Compute negative and positive mass for rho_v, and set negative masses to zero
    if (rho_v(k,j,i,iens) < 0) {
      atomicAdd( rho_v_neg_mass(j,i,iens) , -rho_v(k,j,i,iens) );
      rho_v(k,j,i,iens) = 0;
    } else if (rho_v(k,j,i,iens) > 0) {
      atomicAdd( rho_v_pos_mass(j,i,iens) , rho_v(k,j,i,iens) );
    }
  });

  #ifdef PAM_DEBUG
    ScalarLiveOut<bool> neg_too_large(false);
  #endif

  parallel_for( "Multiplicative hole filler for negative values" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_LAMBDA (int k, int j, int i, int iens) {
    // Redistribute added mass from filling negative values by taking mass away from each positive mass cell
    // Take away mass according to the proportion of total positive mass in this cell
    #ifdef PAM_DEBUG
      if (rho_v_neg_mass(j,i,iens) > rho_v_pos_mass(j,i,iens)) neg_too_large = true;
    #endif
    real factor = rho_v(k,j,i,iens) / rho_v_pos_mass(j,i,iens);
    rho_v(k,j,i,iens) = max( 0._fp , rho_v(k,j,i,iens) - rho_v_neg_mass(j,i,iens) * factor );
  });

  #ifdef PAM_DEBUG
    if (neg_too_large) std::cout << "WARNING: Negative values larger than positive values can fill\n";
  #endif
}


