
#pragma once

#include "pam_const.h"
#include "DataManager.h"

inline void allocate_coupler_state( int nz, int ny, int nx, int nens , DataManager &dm ) {
  dm.register_and_allocate<real>( "density_dry"  , "dry density"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "uvel"         , "x-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "vvel"         , "y-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "wvel"         , "z-direction velocity" , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "temp"         , "temperature"          , {nz,ny,nx,nens} , {"z","y","x","nens"} );
  dm.register_and_allocate<real>( "pressure_dry" , "dry pressure"         , {nz,ny,nx,nens} , {"z","y","x","nens"} );

  auto density_dry  = dm.get_collapsed<real>("density_dry" );
  auto uvel         = dm.get_collapsed<real>("uvel"        );
  auto vvel         = dm.get_collapsed<real>("vvel"        );
  auto wvel         = dm.get_collapsed<real>("wvel"        );
  auto temp         = dm.get_collapsed<real>("temp"        );
  auto pressure_dry = dm.get_collapsed<real>("pressure_dry");

  parallel_for( Bounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
    density_dry (i) = 0;
    uvel        (i) = 0;
    vvel        (i) = 0;
    wvel        (i) = 0;
    temp        (i) = 0;
    pressure_dry(i) = 0;
  });

}


