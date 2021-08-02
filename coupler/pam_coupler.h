
#pragma once

#include "pam_const.h"
#include "DataManager.h"

class PamCoupler {
  public:

  real R_d;
  real R_v;
  real hydrostatic_press_const[4];

  DataManager dm;


  PamCoupler() {
    R_d = 287.;
    R_v = 461.;
    hydrostatic_press_const[0] = -1;
    hydrostatic_press_const[1] = -1;
    hydrostatic_press_const[2] = -1;
    hydrostatic_press_const[3] = -1;
  }


  PamCoupler(real R_d, real R_v) {
    this->R_d = R_d;
    this->R_v = R_v;
    hydrostatic_press_const[0] = -1;
    hydrostatic_press_const[1] = -1;
    hydrostatic_press_const[2] = -1;
    hydrostatic_press_const[3] = -1;
  }


  inline void set_gas_constants(real R_d, real R_v) {
    this->R_d = R_d;
    this->R_v = R_v;
  }


  inline void set_vertical_grid(real2d const &zint_in) {
    int nz   = dm.get_dimension_size("z");
    int nens = dm.get_dimension_size("nens");
    auto zint = dm.get<real,2>("vertical_interface_height");
    auto zmid = dm.get<real,2>("vertical_midpoint_height" );
    parallel_for( SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
      zint(k,iens) = zint_in(k,iens);
      if (k < nz) zmid(k,iens) = 0.5_fp * (zint_in(k,iens) + zint_in(k+1,iens));
    });
  }


  inline void set_vertical_grid(real1d const &zint_in) {
    int nz   = dm.get_dimension_size("z");
    int nens = dm.get_dimension_size("nens");
    auto zint = dm.get<real,2>("vertical_interface_height");
    auto zmid = dm.get<real,2>("vertical_midpoint_height" );
    parallel_for( SimpleBounds<2>(nz+1,nens) , YAKL_LAMBDA (int k, int iens) {
      zint(k,iens) = zint_in(k);
      if (k < nz) zmid(k,iens) = 0.5_fp * (zint_in(k) + zint_in(k+1));
    });
  }


  inline void allocate_coupler_state( int nz, int ny, int nx, int nens ) {
    dm.register_and_allocate<real>( "density_dry"               , "dry density"               , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "uvel"                      , "x-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "vvel"                      , "y-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "wvel"                      , "z-direction velocity"      , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "temp"                      , "temperature"               , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "diag_press"                , "pressure"                  , {nz,ny,nx,nens} , {"z","y","x","nens"} );
    dm.register_and_allocate<real>( "vertical_interface_height" , "vertical interface height" , {nz+1    ,nens} , {"zp1"      ,"nens"} );
    dm.register_and_allocate<real>( "vertical_midpoint_height"  , "vertical midpoint height"  , {nz      ,nens} , {"z"        ,"nens"} );

    auto density_dry  = dm.get_collapsed<real>("density_dry"              );
    auto uvel         = dm.get_collapsed<real>("uvel"                     );
    auto vvel         = dm.get_collapsed<real>("vvel"                     );
    auto wvel         = dm.get_collapsed<real>("wvel"                     );
    auto temp         = dm.get_collapsed<real>("temp"                     );
    auto diag_press   = dm.get_collapsed<real>("diag_press"               );
    auto zint         = dm.get_collapsed<real>("vertical_interface_height");
    auto zmid         = dm.get_collapsed<real>("vertical_midpoint_height" );

    parallel_for( Bounds<1>(nz*ny*nx*nens) , YAKL_LAMBDA (int i) {
      density_dry (i) = 0;
      uvel        (i) = 0;
      vvel        (i) = 0;
      wvel        (i) = 0;
      temp        (i) = 0;
      diag_press  (i) = 0;
      if (i < (nz+1)*nens) zint(i) = 0;
      if (i < (nz  )*nens) zmid(i) = 0;
    });
  }


  YAKL_INLINE real compute_pressure( real rho_d, real rho_v, real T, real R_d, real R_v ) {
    return rho_d*R_d*T + rho_v*R_v*T;
  }


  inline void update_diagnostic_pressure( ) {
    auto dens_dry = dm.get_lev_col<real>("density_dry");
    auto dens_wv  = dm.get_lev_col<real>("water_vapor");
    auto temp     = dm.get_lev_col<real>("temperature");
    auto pressure = dm.get_lev_col<real>("diag_press" );

    int nz   = dens_dry.dimension[0];
    int ncol = dens_dry.dimension[1];

    parallel_for( SimpleBounds<2>(nz,ncol) , YAKL_LAMBDA (int k, int i) {
      real rho_d = dens_dry(k,i);
      real rho_v = dens_wv (k,i);
      real T     = temp    (k,i);
      pressure(k,i) = compute_pressure( rho_d , rho_v , T , R_d , R_v );
    });
  }


  inline void update_hydrostatic_profile( ) {
    auto pressure = dm.get_lev_col<real>("diag_press" );

    int nz   = pressure.dimension[0];
    int ncol = pressure.dimension[1];
  }

};


