
#pragma once

#include "awfl_const.h"
#include "pam_coupler.h"
#include "compute_time_step.h"
#include "reconstruct.h"
#include "tendencies_rho_theta.h"
#include "tendencies_rho_theta_compressible_explicit_semidiscrete.h"

class Dycore_rho_theta_compressible_explicit_ader {
public:
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");

  real dtInit;       // Initial time step (used throughout the simulation)
  awfl::Recon recon;
  awfl::tendencies_rho_theta::Hydrostasis hydrostasis;



  void init(pam::PamCoupler &coupler) {
    if (! coupler.tracer_exists("water_vapor")) endrun("ERROR: processed registered tracers, and water_vapor was not found");
    dtInit = 0;   // Inialize time step to zero, and dimensional splitting switch
    recon      .allocate_and_initialize(coupler);
    hydrostasis.allocate_and_initialize(coupler);
    awfl::tendencies_rho_theta::init_idealized_state_and_tracers( coupler );
  }



  real compute_time_step(pam::PamCoupler const &coupler, real cfl = 0.5) {
    if (dtInit == 0) { dtInit = awfl::compute_time_step(coupler,cfl); }
    return dtInit;
  }



  void timeStep( pam::PamCoupler &coupler , real dtphys ) {
    using awfl::tendencies_rho_theta::compressible_explicit_semidiscrete::compute_tendencies_xyz;
    using awfl::tendencies_rho_theta::compute_mass;
    using yakl::intrinsics::sum;
    using yakl::intrinsics::abs;
    using yakl::intrinsics::maxval;
    using yakl::componentwise::operator*;
    using yakl::componentwise::operator+;
    using yakl::componentwise::operator-;
    using yakl::componentwise::operator/;
    using awfl::tendencies_rho_theta::idR;

    auto state       = awfl::tendencies_rho_theta::createStateArr (coupler);
    auto tracers     = awfl::tendencies_rho_theta::createTracerArr(coupler);
    auto num_state   = state  .extent(0);
    auto num_tracers = tracers.extent(0);
    auto tracer_pos  = coupler.get_tracer_positivity_array();
    auto nx          = coupler.get_nx();
    auto ny          = coupler.get_ny();
    auto nz          = coupler.get_nz();
    auto nens        = coupler.get_nens();

    awfl::tendencies_rho_theta::convert_coupler_state_to_dynamics( coupler , state , tracers , hydrostasis );

    real dt = compute_time_step( coupler );
    int n_iter = ceil( dtphys / dt );
    dt = dtphys / n_iter;

    for (int iter = 0; iter < n_iter; iter++) {

      #ifdef PAM_DEBUG
        auto mass_init = compute_mass( coupler , state , tracers );
      #endif
      auto state_tend  = awfl::tendencies_rho_theta::createStateTendArr (coupler);
      auto tracer_tend = awfl::tendencies_rho_theta::createTracerTendArr(coupler);
      auto state_tmp   = state  .createDeviceObject();
      auto tracers_tmp = tracers.createDeviceObject();

      ///////////////////
      // Stage 1
      ///////////////////
      compute_tendencies_xyz( coupler , state , state_tend , tracers , tracer_tend , recon , hydrostasis , dt );
      parallel_for( "Temporal_ader.h apply tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          state_tmp(l,hs+k,hs+j,hs+i,iens) = state(l,hs+k,hs+j,hs+i,iens) + dt * state_tend(l,k,j,i,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers_tmp(l,hs+k,hs+j,hs+i,iens) = tracers(l,hs+k,hs+j,hs+i,iens) + dt * tracer_tend(l,k,j,i,iens);
          if (tracer_pos(l)) tracers_tmp(l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers_tmp(l,hs+k,hs+j,hs+i,iens) );
        }
      });

      ///////////////////
      // Stage 2
      ///////////////////
      compute_tendencies_xyz( coupler , state_tmp , state_tend , tracers_tmp , tracer_tend , recon , hydrostasis , dt );
      parallel_for( "Temporal_ader.h apply tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          state_tmp(l,hs+k,hs+j,hs+i,iens) = 3._fp*state     (l,hs+k,hs+j,hs+i,iens)/4._fp +
                                                   state_tmp (l,hs+k,hs+j,hs+i,iens)/4._fp + 
                                                dt*state_tend(l,   k,   j,   i,iens)/4._fp;
        }
        for (int l=0; l < num_tracers; l++) {
          tracers_tmp(l,hs+k,hs+j,hs+i,iens) = 3._fp*tracers    (l,hs+k,hs+j,hs+i,iens)/4._fp +
                                                     tracers_tmp(l,hs+k,hs+j,hs+i,iens)/4._fp + 
                                                  dt*tracer_tend(l,   k,   j,   i,iens)/4._fp;
          if (tracer_pos(l)) tracers_tmp(l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers_tmp(l,hs+k,hs+j,hs+i,iens) );
        }
      });

      ///////////////////
      // Stage 3
      ///////////////////
      compute_tendencies_xyz( coupler , state_tmp , state_tend , tracers_tmp , tracer_tend , recon , hydrostasis , dt );
      parallel_for( "Temporal_ader.h apply tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                    YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          state(l,hs+k,hs+j,hs+i,iens) =       state     (l,hs+k,hs+j,hs+i,iens)/3._fp +
                                         2._fp*state_tmp (l,hs+k,hs+j,hs+i,iens)/3._fp + 
                                      2._fp*dt*state_tend(l,   k,   j,   i,iens)/3._fp;
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,hs+i,iens) =       tracers    (l,hs+k,hs+j,hs+i,iens)/3._fp +
                                           2._fp*tracers_tmp(l,hs+k,hs+j,hs+i,iens)/3._fp + 
                                        2._fp*dt*tracer_tend(l,   k,   j,   i,iens)/3._fp;
          if (tracer_pos(l)) tracers(l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers(l,hs+k,hs+j,hs+i,iens) );
        }
      });

      #ifdef PAM_DEBUG
        {
          auto mass_rel_diff = abs(compute_mass( coupler , state , tracers ) - mass_init) / (mass_init + 1.e-20);
          if (maxval(mass_rel_diff) > 1.e-12) { std::cout << mass_rel_diff; endrun("ERROR: Mass not conserved"); }
        }
      #endif

      if (coupler.get_num_dycore_functions() > 0) {
        awfl::tendencies_rho_theta::convert_dynamics_to_coupler_state( coupler , state , tracers );
        coupler.run_dycore_functions( dt );
        awfl::tendencies_rho_theta::convert_coupler_state_to_dynamics( coupler , state , tracers , hydrostasis );
      }

    }

    awfl::tendencies_rho_theta::convert_dynamics_to_coupler_state( coupler , state , tracers );
  }



  void finalize(pam::PamCoupler &coupler) { }



  const char * dycore_name() const { return "AWFL"; }

};


