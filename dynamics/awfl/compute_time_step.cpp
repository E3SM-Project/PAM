
#include "compute_time_step.h"

namespace awfl {

  // Given the model data and CFL value, compute the maximum stable time step
  real compute_time_step(pam::PamCoupler const &coupler, real cfl) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    // TODO: Grab coupler attribute for "anelastic" or "hevi" for alternative calculations

    auto &dm = coupler.get_data_manager_readonly();

    auto dx    = coupler.get_dx();
    auto dy    = coupler.get_dy();
    auto dz    = dm.get<real const,2>("vertical_cell_dz");
    auto nx    = coupler.get_nx();
    auto ny    = coupler.get_ny();
    auto nz    = coupler.get_nz();
    auto nens  = coupler.get_nens();
    auto Rd    = coupler.get_R_d();
    auto cp_d  = coupler.get_cp_d();
    auto cv_d  = cp_d - Rd;
    auto gamma = cp_d / cv_d;
    auto Rv    = coupler.get_R_v();
    auto sim2d = ny == 1;

    // Convert data from DataManager to state and tracers array for convenience
    auto dm_dens_dry = dm.get<real const,4>( "density_dry" );
    auto dm_uvel     = dm.get<real const,4>( "uvel"        );
    auto dm_vvel     = dm.get<real const,4>( "vvel"        );
    auto dm_wvel     = dm.get<real const,4>( "wvel"        );
    auto dm_temp     = dm.get<real const,4>( "temp"        );
    auto dm_dens_vap = dm.get<real const,4>( "water_vapor" );

    // Allocate a 3-D array for the max stable time steps (we'll use this for a reduction later)
    real4d dt3d("dt3d",nz,ny,nx,nens);

    // Loop through the cells, calculate the max stable time step for each cell
    parallel_for( "Spatial.h compute_time_step" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      real rho_d = dm_dens_dry(k,j,i,iens);
      real u     = dm_uvel    (k,j,i,iens);
      real v     = dm_vvel    (k,j,i,iens);
      real w     = dm_wvel    (k,j,i,iens);
      real temp  = dm_temp    (k,j,i,iens);
      real rho_v = dm_dens_vap(k,j,i,iens);
      real p = Rd * rho_d * temp + Rv * rho_v * temp;
      // This neglects non-wv mass-adding tracers, but these are small, and their lack only increases cs
      // Thus the resulting time step is conservative w/r to these missing masses, which is more stable
      real cs = sqrt(gamma*p/(rho_v+rho_d));

      // Compute the maximum stable time step in each direction
      real udt = cfl * dx         / std::max( abs(u-cs) , abs(u+cs) );
      real vdt = cfl * dy         / std::max( abs(v-cs) , abs(v+cs) );
      if (sim2d) vdt = std::numeric_limits<real>::max();
      real wdt = cfl * dz(k,iens) / std::max( abs(w-cs) , abs(w+cs) );

      // Compute the min of the max stable time steps
      dt3d(k,j,i,iens) = std::min( std::min(udt,vdt) , wdt );
    });
    return yakl::intrinsics::minval( dt3d );
  }

}

