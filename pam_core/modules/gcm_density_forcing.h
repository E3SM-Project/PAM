
#pragma once

#include "pam_coupler.h"


namespace modules {

  // Force column-averaged CRM density very strictly to match the GCM density
  inline void gcm_density_forcing( PamCoupler &coupler , real dt ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;

    auto &dm = coupler.get_data_manager_readwrite();

    int nz   = dm.get_dimension_size("z"   );
    int ny   = dm.get_dimension_size("y"   );
    int nx   = dm.get_dimension_size("x"   );
    int nens = dm.get_dimension_size("nens");

    // Current CRM state
    auto gcm_rho_d = dm.get<real,2>( "gcm_density_dry" );
    auto rho_d     = dm.get<real,4>( "density_dry"     );

    // Compute column-average density for current CRM state
    real2d rho_d_colavg("rho_d_colavg",nz,nens);

    memset( rho_d_colavg , 0._fp );

    real r_nx_ny = 1._fp / (nx*ny);  // Avoid costly divisions
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      atomicAdd( rho_d_colavg(k,iens) , rho_d(k,j,i,iens)*r_nx_ny );
    });

    // Apply the difference between GCM and average column densities as forcing
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      // rho_d(k,j,i,iens) /= rho_d_colavg(k,iens);
      // rho_d(k,j,i,iens) *= gcm_rho_d   (k,iens);
      rho_d(k,j,i,iens) += gcm_rho_d(k,iens) - rho_d_colavg(k,iens);
    });
  }

}


