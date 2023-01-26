#pragma once

#include "pam_coupler.h"

class Radiation {
public:

  // Set constants if needed
  Radiation() {
  }

  std::string radiation_name() const {
    return "forced";
  }

  void init(pam::PamCoupler &coupler) {
    auto &dm = coupler.get_data_manager_device_readwrite();
    coupler.set_option<std::string>("radiation","forced");
    auto nens   = coupler.get_option<int>("ncrms");
    auto nz     = coupler.get_option<int>("crm_nz");
    auto rad_nx = coupler.get_option<int>("rad_nx");
    auto rad_ny = coupler.get_option<int>("rad_ny");
    dm.register_and_allocate<real>("rad_enthalpy_tend" ,"radiation tendency from external calculation",{nz,rad_ny,rad_nx,nens},{"z","rad_y","rad_x","nens"});
  }

  void timeStep( pam::PamCoupler &coupler ) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;
    auto &dm = coupler.get_data_manager_device_readwrite();
    auto nens   = coupler.get_option<int>("ncrms");
    auto nz     = coupler.get_option<int>("crm_nz");
    auto crm_nx = coupler.get_option<int>("crm_nx");
    auto crm_ny = coupler.get_option<int>("crm_ny");
    auto rad_nx = coupler.get_option<int>("rad_nx");
    auto rad_ny = coupler.get_option<int>("rad_ny");
    auto dt     = coupler.get_option<real>("crm_dt");
    auto cp_d   = coupler.get_option<real>("cp_d");
    auto rad_enthalpy_tend = dm.get<real const,4>("rad_enthalpy_tend");
    auto temperature       = dm.get<real,4>("temp");
    parallel_for( YAKL_AUTO_LABEL() , SimpleBounds<4>(nz,crm_ny,crm_nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      int i_rad = i / (crm_nx/rad_nx);
      int j_rad = j / (crm_ny/rad_ny);
      temperature(k,j,i,iens) += rad_enthalpy_tend(k,j_rad,i_rad,iens) / cp_d * dt;
    });
  }

  void finalize(pam::PamCoupler &coupler) {
  }

};

