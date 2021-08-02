
#pragma once

#include "pam_const.h"
#include "DataManager.h"

inline void allocate_coupler_state( int nz, int ny, int nx, int nens , DataManager &dm ) {
  dm.register_and_allocate<real>( "density_dry"  , "dry density"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "uvel"         , "x-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "vvel"         , "y-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "wvel"         , "z-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "temp"         , "temperature"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );

  auto density_dry  = dm.get_collapsed<real>("density_dry" );
  auto uvel         = dm.get_collapsed<real>("uvel"        );
  auto vvel         = dm.get_collapsed<real>("vvel"        );
  auto wvel         = dm.get_collapsed<real>("wvel"        );
  auto temp         = dm.get_collapsed<real>("temp"        );

  parallel_for( Bounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
    density_dry (i) = 0;
    uvel        (i) = 0;
    vvel        (i) = 0;
    wvel        (i) = 0;
    temp        (i) = 0;
  });
}


YAKL_INLINE real compute_pressure( real rho_d, real rho_v, real T, real R_d, real R_v ) {
  return rho_d*R_d*T + rho_v*R_v*T;
}


YAKL_INLINE real4d compute_pressure( real4d const &dens_dry, real4d const &dens_wv, real4d const &temp,
                                     real R_d, real R_v ) {
  int nz   = dens_dry.dimension[0];
  int ny   = dens_dry.dimension[1];
  int nx   = dens_dry.dimension[2];
  int nens = dens_dry.dimension[3];

  real4d pressure("pressure",nz,ny,nx,nens);

  parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
    real rho_d = dens_dry(k,j,i,iens);
    real rho_v = dens_wv (k,j,i,iens);
    real T     = temp    (k,j,i,iens);
    pressure(k,j,i,iens) = compute_pressure( rho_d , rho_v , T , R_d , R_v );
  });

  return pressure;
}


