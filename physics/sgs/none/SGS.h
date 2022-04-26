
#pragma once


int static constexpr num_tracers_sgs = 0;

#include "pam_coupler.h"


class SGS {
public:

  // TODO: Make sure the constants vibe with P3
  // Set constants and likely num_tracers as well, and anything else you can do immediately
  SGS() {
  }



  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(pam::PamCoupler &coupler) {
    coupler.set_option<std::string>("sgs","none");
  }



  void timeStep( pam::PamCoupler &coupler , real dt ) {
  }



  std::string sgs_name() const {
    return "none";
  }


  void finalize(pam::PamCoupler &coupler) {
  }



};



