
#pragma once

#include "pam_coupler.h"

class Microphysics {
public:

  real R_d    ;
  real cp_d   ;
  real cv_d   ;
  real gamma_d;
  real kappa_d;
  real R_v    ;
  real cp_v   ;
  real cv_v   ;
  real p0     ;
  real grav   ;



  Microphysics() {
    R_d     = 287.;
    cp_d    = 1003.;
    cv_d    = cp_d - R_d;
    gamma_d = cp_d / cv_d;
    kappa_d = R_d  / cp_d;
    R_v     = 461.;
    cp_v    = 1859;
    cv_v    = R_v - cp_v;
    p0      = 1.e5;
    grav    = 9.81;
  }



  // Gotta have at least three tracers- the mass of water species
  static int constexpr get_num_tracers() {
    return 1;
  }



  // Have to declare at least water vapor
  void init(pam::PamCoupler &coupler) {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    int nx   = coupler.get_nx  ();
    int ny   = coupler.get_ny  ();
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    //                 name              description            positive   adds mass
    coupler.add_tracer("water_vapor"  , "Water Vapor"   ,       true     , true);

    auto &dm = coupler.get_data_manager_device_readwrite();

    // Zero out the tracers

    auto rho_v = dm.get_collapsed<real>("water_vapor");
    parallel_for("zero water vapor", nz*ny*nx*nens , YAKL_LAMBDA (int i) { rho_v(i) = 0; } );

    coupler.set_option<std::string>("micro","none");
    coupler.set_option<real>("R_d" ,R_d  );
    coupler.set_option<real>("R_v" ,R_v  );
    coupler.set_option<real>("cp_d",cp_d );
    coupler.set_option<real>("cp_v",cp_v );
    coupler.set_option<real>("grav",grav );
    coupler.set_option<real>("p0"  ,p0   );
  }



  void timeStep( pam::PamCoupler &coupler ) {
    // Do microphysicsy stuff to the coupler variables and the tracers
  }



  std::string micro_name() const {
    return "none";
  }

  void finalize(pam::PamCoupler &coupler) {
  }
};
