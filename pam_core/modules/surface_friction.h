
#pragma once

#include "pam_coupler.h"

namespace modules {

  real constexpr vonk = 0.4;
  real constexpr eps = 1.0e-10;
  real constexpr am = 4.8;
  real constexpr bm = 19.3;
  real constexpr pi = 3.14159;

  // estimate roughness height for momentum
  YAKL_INLINE void z0_est(real z, real bflx, real wnd, real ustar, real &z0) {
    real c1 = pi/2.0 - 3.0*log(2.0);
    real rlmo = -bflx*vonk/(ustar*ustar*ustar+eps);
    real zeta = std::min(1.0,z*rlmo);
    real x;
    real psi1;
    if(zeta >= 0.0) {
      psi1 = -am*zeta;
    }
    else {
      x = sqrt(sqrt(1.0-bm*zeta));
      psi1 = 2.0*log(1.0+x) + log(1.0+x*x) -2.0*atan(x) + c1;
    }
    real lnz = std::max(0.0, vonk*wnd/(ustar+eps) +psi1);
    z0 = z*exp(-lnz);
  }


  // calculate friction speed (sqrt of momentum flux) [m/s]
  YAKL_INLINE real diag_ustar(real z, real bflx, real wnd, real z0) {
    real lnz = log(z/z0);
    real klnz = vonk/lnz;
    real c1 = pi/2.0 - 3.0*log(2.0);
    real ustar = wnd*klnz;
    if (bflx != 0.0) {
      for (int iterate = 0; iterate < 8; iterate++) {
        real rlmo = -bflx * vonk/(ustar*ustar*ustar + eps);
        real zeta = std::min(1.0,z*rlmo);
        if (zeta>0.0) {
          ustar = vonk*wnd / (lnz+am*zeta);
        } else {
          real x = sqrt(sqrt(1.0-bm*zeta));
          real psi1 = 2.0*log(1.0+x) + log(1.0+x*x) - 2.0*atan(x) + c1;
          ustar = wnd*vonk/(lnz-psi1);
        }
      }
    }
    return ustar;
  }


  inline void surface_friction_init( pam::PamCoupler &coupler, realConst1d &tau_in, realConst1d &bflx_in ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;
    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();
    auto &dm = coupler.get_data_manager_device_readwrite();
    dm.register_and_allocate<real>( "z0"      , "Momentum roughness height [m]",     {nens},{"nens"} );
    dm.register_and_allocate<real>( "sfc_bflx", "large-scale sfc buoyancy flux [K m/s]", {nens},{"nens"} );
    auto z0       = dm.get<real,1>("z0");
    auto rho_d    = dm.get<real,4>("density_dry");
    auto rho_v    = dm.get<real,4>("water_vapor");
    auto zmid     = dm.get<real,2>("vertical_midpoint_height" );
    auto sfc_bflx = dm.get<real,1>("sfc_bflx");
    auto gcm_uvel = dm.get<real,2>("gcm_uvel");
    auto gcm_vvel = dm.get<real,2>("gcm_vvel");

    real1d rho_horz_mean("rho_horz_mean",nens);
    real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
    parallel_for("compute horz mean wind", SimpleBounds<3>(ny,nx,nens), YAKL_LAMBDA (int j, int i, int iens) {
      atomicAdd( rho_horz_mean(iens), ( rho_d(0,j,i,iens) + rho_v(0,j,i,iens) ) * r_nx_ny );
    });

    parallel_for( nens , YAKL_LAMBDA (int iens) {
      sfc_bflx(iens) = bflx_in(iens);
      real wnd_spd = std::max( 1.0, sqrt( gcm_uvel(0,iens)*gcm_uvel(0,iens) 
                                    +gcm_vvel(0,iens)*gcm_vvel(0,iens) ) );
      real ustar = sqrt( tau_in(iens) / rho_horz_mean(iens) );
      z0_est( zmid(0,iens), sfc_bflx(iens), wnd_spd, ustar, z0(iens));
      z0(iens) = std::max(0.00001,std::min(1.0,z0(iens)));
    });
  }


  inline void compute_surface_friction( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    using yakl::atomicAdd;
    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();
    auto &dm = coupler.get_data_manager_device_readwrite();
    
    auto z0            = dm.get<real,1>("z0");             // momentum roughness height
    auto sfc_bflx      = dm.get<real,1>("sfc_bflx");       // surface buoyancy flux
    auto sfc_mom_flx_u = dm.get<real,3>("sfc_mom_flx_u" ); // momentum fluxes applied in SGS scheme
    auto sfc_mom_flx_v = dm.get<real,3>("sfc_mom_flx_v" ); // momentum fluxes applied in SGS scheme
    auto zmid          = dm.get<real const,2>("vertical_midpoint_height" );
    auto zint          = dm.get<real const,2>("vertical_interface_height");
    auto rho_d         = dm.get<real const,4>("density_dry");
    auto rho_v         = dm.get<real const,4>("water_vapor");
    auto uvel          = dm.get<real const,4>("uvel");
    auto vvel          = dm.get<real const,4>("vvel");

    real1d u_horz_mean("u_horz_mean",nens);
    real1d v_horz_mean("v_horz_mean",nens);
    real1d rho_horz_mean("rho_horz_mean",nens);
    real r_nx_ny  = 1._fp / (nx*ny);  // precompute reciprocal to avoid costly divisions
    parallel_for("compute horz mean wind", SimpleBounds<3>(ny,nx,nens), YAKL_LAMBDA (int j, int i, int iens) {
      atomicAdd( u_horz_mean(iens), uvel(0,j,i,iens) * r_nx_ny );
      atomicAdd( v_horz_mean(iens), vvel(0,j,i,iens) * r_nx_ny );
      atomicAdd( rho_horz_mean(iens), ( rho_d(0,j,i,iens) + rho_v(0,j,i,iens) ) * r_nx_ny );
    });

    parallel_for( "surface momentum flux" , SimpleBounds<3>(ny,nx,nens) , YAKL_LAMBDA (int j, int i, int iens) {
      // calculate surface momentum flux identical to SAM with units = [kg m/s2]
      real wnd_spd = std::max( 1.0, sqrt( uvel(0,j,i,iens)*uvel(0,j,i,iens) + vvel(0,j,i,iens)*vvel(0,j,i,iens) ) );
      real ustar = diag_ustar( zmid(0,iens), sfc_bflx(iens), wnd_spd, z0(iens)); 
      real tau00 = rho_horz_mean(iens) * ustar * ustar;
      sfc_mom_flx_u(j,i,iens) = -( uvel(0,j,i,iens) - u_horz_mean(iens) ) / wnd_spd * tau00;
      sfc_mom_flx_v(j,i,iens) = -( vvel(0,j,i,iens) - v_horz_mean(iens) ) / wnd_spd * tau00;
      // convert units to [m2/s2] for SHOC - interpolate to get density at surface
      real rho0 = rho_d(0,j,i,iens) + rho_v(0,j,i,iens);
      real rho1 = rho_d(1,j,i,iens) + rho_v(1,j,i,iens);
      real rho_sfc = 2.0*rho0 - rho1;
      real dz = zint(1,iens) - zint(0,iens);
      sfc_mom_flx_u(j,i,iens) = sfc_mom_flx_u(j,i,iens) * rho_sfc / dz ;
      sfc_mom_flx_v(j,i,iens) = sfc_mom_flx_v(j,i,iens) * rho_sfc / dz ;
    });
  }


}


