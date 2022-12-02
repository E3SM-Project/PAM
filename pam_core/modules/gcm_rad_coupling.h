
#pragma once

#include "pam_coupler.h"

namespace modules {
  
  void apply_rad_tendencies( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    
    int nz   = coupler.get_nz  ();
    int ny   = coupler.get_ny  ();
    int nx   = coupler.get_nx  ();
    int nens = coupler.get_nens();
    
    auto &dm = coupler.get_data_manager_readwrite();

    auto dt = coupler.get_option<real>("crm_dt");

    auto nx_rad = coupler.get_option<int>("nx_rad");
    auto ny_rad = coupler.get_option<int>("ny_rad");
    auto cp_d   = coupler.get_cp_d();
    auto rad_enthalpy_tend = dm.get<real const,4>("crm_rad_tend");
    auto temp = dm.get<real,4>("temp");

    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      int i_rad = i / (nx/nx_rad);
      int j_rad = j / (nx/ny_rad);
      temp(k,j,i,iens) += rad_enthalpy_tend(k,j_rad,i_rad,iens) / cp_d * dt;
    });
  }

}


