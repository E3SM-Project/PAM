
#pragma once

//#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"

using pam::PamCoupler;

int static constexpr num_tracers_micro = 1;

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
  int get_num_tracers() const {
    return 1;
  }



  // Have to declare at least water species
  void init(PamCoupler &coupler) {
    int nx   = coupler.get_nx  ();
    int ny   = coupler.get_ny  ();
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    //                 name              description            positive   adds mass
    coupler.add_tracer("water_vapor"  , "Water Vapor"   ,       true     , true);

    // Zero out the tracers
    auto rho_v = coupler.dm.get_collapsed<real>("water_vapor");

    parallel_for( nz*ny*nx*nens , YAKL_LAMBDA (int i) { 
	rho_v(i) = 0; 
	} );

    coupler.set_option<std::string>("micro","none");
  }



  void timeStep( PamCoupler &coupler , real dt ) {
    // Do microphysicsy stuff to the coupler variables and the tracers
  }



  std::string micro_name() const {
    return "none";
  }


};
