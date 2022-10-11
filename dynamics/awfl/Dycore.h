
#pragma once

#include "awfl_const.h"
#include "pam_coupler.h"
#include "compute_time_step.h"
#include "reconstruct.h"
#include "tendencies_rho_theta.h"
#include "tendencies_rho_theta_compressible_explicit_ader.h"

class Dycore {
public:
  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");

  bool dim_switch;
  real dtInit;       // Initial time step (used throughout the simulation)
  awfl::Recon recon;
  awfl::tendencies_rho_theta::Hydrostasis hydrostasis;



  void init(pam::PamCoupler &coupler) {
    if (! coupler.tracer_exists("water_vapor")) endrun("ERROR: processed registered tracers, and water_vapor was not found");
    dtInit = 0;   // Inialize time step to zero, and dimensional splitting switch
    recon      .allocate_and_initialize(coupler);
    hydrostasis.allocate_and_initialize(coupler);
    awfl::tendencies_rho_theta::init_idealized_state_and_tracers( coupler );
    dim_switch = true;
  }



  real compute_time_step(pam::PamCoupler const &coupler, real cfl = 0.75) {
    if (dtInit == 0) { dtInit = awfl::compute_time_step(coupler,cfl); }
    return dtInit;
  }



  void apply_tendencies(real5d const &state   , realConst5d state_tend  ,
                        real5d const &tracers , realConst5d tracer_tend , real dt ) {
    int num_state   = state_tend.extent(0);
    int nz          = state_tend.extent(1);
    int ny          = state_tend.extent(2);
    int nx          = state_tend.extent(3);
    int nens        = state_tend.extent(4);
    int num_tracers = tracer_tend.extent(0);

    parallel_for( "Temporal_ader.h apply tendencies" , SimpleBounds<4>(nz,ny,nx,nens) ,
                  YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        state(l,hs+k,hs+j,hs+i,iens) += dt * state_tend(l,k,j,i,iens);
      }
      for (int l=0; l < num_tracers; l++) {
        tracers(l,hs+k,hs+j,hs+i,iens) += dt * tracer_tend(l,k,j,i,iens);
        tracers(l,hs+k,hs+j,hs+i,iens) = std::max( 0._fp , tracers(l,hs+k,hs+j,hs+i,iens) );
      }
    });
  }



  void timeStep( pam::PamCoupler &coupler , real dtphys ) {
    using awfl::tendencies_rho_theta::compressible_explicit_ader::compute_tendencies_x;
    using awfl::tendencies_rho_theta::compressible_explicit_ader::compute_tendencies_y;
    using awfl::tendencies_rho_theta::compressible_explicit_ader::compute_tendencies_z;

    auto state       = awfl::tendencies_rho_theta::createStateArr     (coupler);
    auto tracers     = awfl::tendencies_rho_theta::createTracerArr    (coupler);

    awfl::tendencies_rho_theta::convert_coupler_state_to_dynamics( coupler , state , tracers , hydrostasis );

    real dt = compute_time_step( coupler );
    int n_iter = ceil( dtphys / dt );
    dt = dtphys / n_iter;

    for (int iter = 0; iter < n_iter; iter++) {

      #ifdef PAM_DEBUG
        validate_array_positive(tracers);
        validate_array_inf_nan(state);
        validate_array_inf_nan(tracers);
        std::vector<real> mass_init = compute_mass( coupler , state , tracers );
      #endif
      auto state_tend  = awfl::tendencies_rho_theta::createStateTendArr (coupler);
      auto tracer_tend = awfl::tendencies_rho_theta::createTracerTendArr(coupler);

      if (dim_switch) {
        compute_tendencies_x( coupler , state , state_tend , tracers , tracer_tend , recon , dt );
        apply_tendencies    ( state , state_tend , tracers , tracer_tend , dt );

        compute_tendencies_y( coupler , state , state_tend , tracers , tracer_tend , recon , dt );
        apply_tendencies    ( state , state_tend , tracers , tracer_tend , dt );

        compute_tendencies_z( coupler , state , state_tend , tracers , tracer_tend , recon , hydrostasis , dt );
        apply_tendencies    ( state , state_tend , tracers , tracer_tend , dt );
      } else {
        compute_tendencies_z( coupler , state , state_tend , tracers , tracer_tend , recon , hydrostasis , dt );
        apply_tendencies    ( state , state_tend , tracers , tracer_tend , dt );

        compute_tendencies_y( coupler , state , state_tend , tracers , tracer_tend , recon , dt );
        apply_tendencies    ( state , state_tend , tracers , tracer_tend , dt );

        compute_tendencies_x( coupler , state , state_tend , tracers , tracer_tend , recon , dt );
        apply_tendencies    ( state , state_tend , tracers , tracer_tend , dt );
      }
      dim_switch = ! dim_switch;

      #ifdef PAM_DEBUG
        std::vector<real> mass_final = compute_mass( coupler , state , tracers );
        for (int l=0; l < mass_final.size(); l++) {
          real mass_diff;
          if (mass_init[l] > 0) { mass_diff = abs(mass_final[l] - mass_init[l]) / abs(mass_init[l]); }
          else                  { mass_diff = mass_final[l]; }
          real tol = std::is_same<real,float>::value ? 1.e-5 : 1.e-12;
          if (mass_diff > tol) { endrun("ERROR: mass not conserved by dycore"); }
        }
        validate_array_positive(tracers);
        validate_array_inf_nan(state);
        validate_array_inf_nan(tracers);
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


