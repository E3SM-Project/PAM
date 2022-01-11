
#include "pam_const.h"
#include "pam_coupler.h"


inline void compute_gcm_forcing_tendencies( PamCoupler &coupler , real2d &rho_d_in , real2d &uvel_in , real2d &vvel_in ,
                                            real2d &wvel_in , real2d &temp_in , real2d &rho_v_in , real dt_gcm ) {
  using yakl::atomicAdd;
  auto &dm = coupler.dm;
  // Get current state from coupler
  auto rho_d = dm.get<real,4>( "density_dry" );
  auto uvel  = dm.get<real,4>( "uvel"        );
  auto vvel  = dm.get<real,4>( "vvel"        );
  auto wvel  = dm.get<real,4>( "wvel"        );
  auto temp  = dm.get<real,4>( "temp"        );
  auto rho_v = dm.get<real,4>( "water_vapor" );

  int nz   = dm.get_dimension("z"   );
  int ny   = dm.get_dimension("y"   );
  int nx   = dm.get_dimension("x"   );
  int nens = dm.get_dimension("nens");

  real2d colsum_rho_d("colsum_rho_d",nz,nens);
  real2d colsum_uvel ("colsum_uvel" ,nz,nens);
  real2d colsum_vvel ("colsum_vvel" ,nz,nens);
  real2d colsum_wvel ("colsum_wvel" ,nz,nens);
  real2d colsum_temp ("colsum_temp" ,nz,nens);
  real2d colsum_rho_v("colsum_rho_v",nz,nens);

  parallel_for( "Initialize colsum to zero" , SimpleBounds<2>(nz,nens) , 
                YAKL_LAMBDA (int k, int iens) {
    colsum_rho_d(k,iens) = 0;
    colsum_uvel (k,iens) = 0;
    colsum_vvel (k,iens) = 0;
    colsum_wvel (k,iens) = 0;
    colsum_temp (k,iens) = 0;
    colsum_rho_v(k,iens) = 0;
  });

  parallel_for( "Compute summed column of current state" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_LAMBDA (int k, int j, int i, int iens) {
    atomicAdd( colsum_rho_d(k,iens) , rho_d(k,j,i,iens) );
    atomicAdd( colsum_uvel (k,iens) , uvel (k,j,i,iens) );
    atomicAdd( colsum_vvel (k,iens) , vvel (k,j,i,iens) );
    atomicAdd( colsum_wvel (k,iens) , wvel (k,j,i,iens) );
    atomicAdd( colsum_temp (k,iens) , temp (k,j,i,iens) );
    atomicAdd( colsum_rho_v(k,iens) , rho_v(k,j,i,iens) );
  });

  // We need the GCM forcing tendencies later, so store these in the coupler's data manager
  dm.register_and_allocate<real>( "gcm_tend_rho_d" , "GCM forcing tendency for dry density"         , {nz,nens} , {"z","nens"} );
  dm.register_and_allocate<real>( "gcm_tend_uvel"  , "GCM forcing tendency for u-velocity"          , {nz,nens} , {"z","nens"} );
  dm.register_and_allocate<real>( "gcm_tend_vvel"  , "GCM forcing tendency for v-velocity"          , {nz,nens} , {"z","nens"} );
  dm.register_and_allocate<real>( "gcm_tend_wvel"  , "GCM forcing tendency for w-velocity"          , {nz,nens} , {"z","nens"} );
  dm.register_and_allocate<real>( "gcm_tend_temp"  , "GCM forcing tendency for temperature"         , {nz,nens} , {"z","nens"} );
  dm.register_and_allocate<real>( "gcm_tend_rho_v" , "GCM forcing tendency for water vapor density" , {nz,nens} , {"z","nens"} );

  // Retrieve the GCM forcing tendency arrays
  auto gcm_tend_rho_d = dm.get<real,2>("gcm_tend_rho_d");
  auto gcm_tend_uvel  = dm.get<real,2>("gcm_tend_uvel" );
  auto gcm_tend_vvel  = dm.get<real,2>("gcm_tend_vvel" );
  auto gcm_tend_wvel  = dm.get<real,2>("gcm_tend_wvel" );
  auto gcm_tend_temp  = dm.get<real,2>("gcm_tend_temp" );
  auto gcm_tend_rho_v = dm.get<real,2>("gcm_tend_rho_v");

  real r_nx_ny  = 1._fp / (nx*ny);
  real r_dt_gcm = 1._fp / dt_gcm;
  parallel_for( "Compute GCM forcing tendencies" , SimpleBounds<2>(nz,nens) , 
                YAKL_LAMBDA (int k, int iens) {
    gcm_tend_rho_d(k,iens) = ( rho_d_in(k,iens) - colsum_rho_d(k,iens)*r_nx_ny ) * r_dt_gcm;
    gcm_tend_uvel (k,iens) = ( uvel_in (k,iens) - colsum_uvel (k,iens)*r_nx_ny ) * r_dt_gcm;
    gcm_tend_vvel (k,iens) = ( vvel_in (k,iens) - colsum_vvel (k,iens)*r_nx_ny ) * r_dt_gcm;
    gcm_tend_wvel (k,iens) = ( wvel_in (k,iens) - colsum_wvel (k,iens)*r_nx_ny ) * r_dt_gcm;
    gcm_tend_temp (k,iens) = ( temp_in (k,iens) - colsum_temp (k,iens)*r_nx_ny ) * r_dt_gcm;
    gcm_tend_rho_v(k,iens) = ( rho_v_in(k,iens) - colsum_rho_v(k,iens)*r_nx_ny ) * r_dt_gcm;
  });

}




inline void apply_gcm_forcing_tendencies( PamCoupler &coupler , real dt ) {
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

  parallel_for( "Compute summed column of current state" , SimpleBounds<4>(nz,ny,nx,nens) , 
                YAKL_LAMBDA (int k, int j, int i, int iens) {
    rho_d(k,j,i,iens) += gcm_tend_rho_d(k,iens) * dt;
    uvel (k,j,i,iens) += gcm_tend_uvel (k,iens) * dt;
    vvel (k,j,i,iens) += gcm_tend_vvel (k,iens) * dt;
    wvel (k,j,i,iens) += gcm_tend_wvel (k,iens) * dt;
    temp (k,j,i,iens) += gcm_tend_temp (k,iens) * dt;
    rho_v(k,j,i,iens) += gcm_tend_rho_v(k,iens) * dt;
  });
}


