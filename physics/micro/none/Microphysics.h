
#pragma once

#include "awfl_const.h"
#include "DataManager.h"
#include "pam_coupler.h"

using pam::PamCoupler;

class Microphysics {
public:

  int tracer_index_vapor;

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



  // Gotta have at least one tracer
  int get_num_tracers() const {
    return 1;
  }



  YAKL_INLINE int get_water_vapor_index() const {
    return tracer_index_vapor;
  }



  // Have to declare at least water vapor
  template <class DC>
  void init(DC &dycore , PamCoupler &coupler) {
    int nx   = coupler.get_nx  ();
    int ny   = coupler.get_ny  ();
    int nz   = coupler.get_nz  ();
    int nens = coupler.get_nens();

    // Register tracers in the dycore
    //                                     name              description       positive   adds mass
    tracer_index_vapor = dycore.add_tracer("water_vapor"   , "Water Vapor"   , true     , true);

    // Register and allocate the tracers in the DataManager
    coupler.dm.register_and_allocate<real>( "water_vapor"   , "Water Vapor"   , {nz,ny,nx,nens} , {"z","y","x","nens"} );

    auto rho_v = coupler.dm.get<real,4>("water_vapor");
    parallel_for( Bounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) { rho_v(k,j,i,iens) = 0; } );
  }



  void timeStep( PamCoupler &coupler , real dt , real etime ) {
    // Do microphysicsy stuff to the coupler variables and the tracers
  }



  std::string micro_name() const {
    return "none";
  }


};
