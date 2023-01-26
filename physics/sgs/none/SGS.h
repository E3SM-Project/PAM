
#pragma once

#include "pam_coupler.h"


class SGS {
public:

  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  SGS() {
  }

  static int constexpr get_num_tracers() {
    // Store surface momentum fluxes in data manager to facilitate internal surface calculations
    dm.register_and_allocate<real>( "sfc_mom_flx_u", "Surface flux of U-momentum", {ny,nx,nens}, {"y","x","nens"} );
    dm.register_and_allocate<real>( "sfc_mom_flx_v", "Surface flux of V-momentum", {ny,nx,nens}, {"y","x","nens"} );
    return 0;
  }

  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(pam::PamCoupler &coupler) {
    coupler.set_option<std::string>("sgs","none");
  }



  void timeStep( pam::PamCoupler &coupler ) {
  }



  std::string sgs_name() const {
    return "none";
  }


  void finalize(pam::PamCoupler &coupler) {
  }



};



