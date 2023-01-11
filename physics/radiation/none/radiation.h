#pragma once

#include "pam_coupler.h"

class Radiation {
public:

  Radiation() {
  }

  std::string radiation_name() const {
    return "none";
  }

  // Can do whatever you want, but mainly for registering tracers and allocating data
  void init(pam::PamCoupler &coupler) {
    coupler.set_option<std::string>("radiation","none");
  }

  void timeStep( pam::PamCoupler &coupler ) {
  }

  void finalize(pam::PamCoupler &coupler) {
  }

};

