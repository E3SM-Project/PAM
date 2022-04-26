
#pragma once

#include "pam_coupler.h"

namespace modules {

  inline void broadcast_initial_gcm_column( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    auto nz   = coupler.get_nz  ();
    auto ny   = coupler.get_ny  ();
    auto nx   = coupler.get_nx  ();
    auto nens = coupler.get_nens();

    auto &dm = coupler.get_data_manager_readwrite();

    auto gcm_rho_d = dm.get<real const,2>("gcm_density_dry");
    auto gcm_uvel  = dm.get<real const,2>("gcm_uvel"       );
    auto gcm_vvel  = dm.get<real const,2>("gcm_vvel"       );
    auto gcm_wvel  = dm.get<real const,2>("gcm_wvel"       );
    auto gcm_temp  = dm.get<real const,2>("gcm_temp"       );
    auto gcm_rho_v = dm.get<real const,2>("gcm_water_vapor");

    auto rho_d = dm.get<real,4>("density_dry");
    auto uvel  = dm.get<real,4>("uvel"       );
    auto vvel  = dm.get<real,4>("vvel"       );
    auto wvel  = dm.get<real,4>("wvel"       );
    auto temp  = dm.get<real,4>("temp"       );
    auto rho_v = dm.get<real,4>("water_vapor");

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      rho_d(k,j,i,iens) = gcm_rho_d(k,iens);
      uvel (k,j,i,iens) = gcm_uvel (k,iens);
      vvel (k,j,i,iens) = gcm_vvel (k,iens);
      wvel (k,j,i,iens) = gcm_wvel (k,iens);
      temp (k,j,i,iens) = gcm_temp (k,iens);
      rho_v(k,j,i,iens) = gcm_rho_v(k,iens);
    });
  }

}


