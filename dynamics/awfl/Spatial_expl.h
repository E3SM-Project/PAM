
#pragma once

#include "awfl_const.h"
#include "pam_coupler.h"
#include "compute_time_step.h"
#include "reconstruct.h"
#include "tendencies_rho_theta.h"


class Spatial_operator {
public:
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");
  real dtInit;       // Initial time step (used throughout the simulation)
  awfl::Recon recon;
  awfl::tendencies_rho_theta::Hydrostasis hydrostasis;

  // When this class is created, initialize num_tracers to zero
  Spatial_operator() { }


  real compute_time_step(pam::PamCoupler const &coupler, real cfl = 0.75) {
    if (dtInit == 0) { dtInit = awfl::compute_time_step(coupler,cfl); }
    return dtInit;
  }


  // Initialize crap needed by recon()
  void init(pam::PamCoupler &coupler) {
    if (! coupler.tracer_exists("water_vapor")) endrun("ERROR: processed registered tracers, and water_vapor was not found");
    dtInit = 0;   // Inialize time step to zero, and dimensional splitting switch
    recon      .allocate_and_initialize(coupler);
    hydrostasis.allocate_and_initialize(coupler);
    awfl::tendencies_rho_theta::init_idealized_state_and_tracers( coupler );
  }


  const char * getName() { return ""; }


  void finalize(real4d const &state , real4d const &tracers) {}

};

