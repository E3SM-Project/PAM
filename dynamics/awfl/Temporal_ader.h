
#pragma once

#include "awfl_const.h"
#include "pam_coupler.h"

template <class Spatial> class Temporal_operator {
public:

  bool dim_switch;
  Spatial space_op;

  void init(pam::PamCoupler &coupler) {
    space_op.init(coupler);
    dim_switch = true;
  }


  void init_idealized_state_and_tracers( pam::PamCoupler &coupler ) {
    space_op.init_idealized_state_and_tracers( coupler );
  }


  real compute_time_step(pam::PamCoupler const &coupler, real cfl_in = 0.75) {
    return space_op.compute_time_step(coupler, cfl_in);
  }


  std::vector<real> compute_mass( pam::PamCoupler const &coupler , realConst5d state , realConst5d tracers ) const {
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    int nz   = coupler.get_nz();
    int ny   = coupler.get_ny();
    int nx   = coupler.get_nx();
    int nens = coupler.get_nens();

    int idR = Spatial::idR;
    int num_tracers = coupler.get_num_tracers();
    auto dz = coupler.get_data_manager_readonly().get<real const,2>("vertical_cell_dz");

    std::vector<real> mass(num_tracers+1);
    real4d tmp("tmp",nz,ny,nx,nens);

    parallel_for( "Temporal_ader.h state mass" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tmp(k,j,i,iens) = state(idR,hs+k,hs+j,hs+i,iens) * dz(k,iens);
    });
    mass[0] = yakl::intrinsics::sum(tmp);

    for (int l=0; l < num_tracers; l++) {
      parallel_for( "Temporal_ader.h tracer mass" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        tmp(k,j,i,iens) = tracers(l,hs+k,hs+j,hs+i,iens) * dz(k,iens);
      });
      mass[l+1] = yakl::intrinsics::sum(tmp);
    }
    return mass;
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
    using yakl::c::parallel_for;
    using yakl::c::SimpleBounds;

    real dt = compute_time_step( coupler );

    auto state       = awfl::tendencies_rho_theta::createStateArr     (coupler);
    auto tracers     = awfl::tendencies_rho_theta::createTracerArr    (coupler);
    auto state_tend  = awfl::tendencies_rho_theta::createStateTendArr (coupler);
    auto tracer_tend = awfl::tendencies_rho_theta::createTracerTendArr(coupler);

    #ifdef PAM_DEBUG
      state   = 0._fp;
      tracers = 0._fp;
    #endif

    awfl::tendencies_rho_theta::convert_coupler_state_to_dynamics( coupler , state , tracers , space_op.hydrostasis );

    int idR         = space_op.idR;
    int idU         = space_op.idU;
    int idV         = space_op.idV;
    int idW         = space_op.idW;
    int idT         = space_op.idT;
    int nx          = coupler.get_nx();
    int ny          = coupler.get_ny();
    int nz          = coupler.get_nz();
    int nens        = coupler.get_nens();

    int n_iter = ceil( dtphys / dt );
    dt = dtphys / n_iter;

    for (int iter = 0; iter < n_iter; iter++) {

      #ifdef PAM_DEBUG
        validate_array_positive(tracers);
        validate_array_inf_nan(state);
        validate_array_inf_nan(tracers);
        std::vector<real> mass_init = compute_mass( coupler , state , tracers );
      #endif

      if (dim_switch) {
        space_op.compute_tendencies_x( coupler , state , state_tend , tracers , tracer_tend , space_op.recon , dt );
        apply_tendencies             ( state , state_tend , tracers , tracer_tend , dt );

        space_op.compute_tendencies_y( coupler , state , state_tend , tracers , tracer_tend , space_op.recon , dt );
        apply_tendencies             ( state , state_tend , tracers , tracer_tend , dt );

        space_op.compute_tendencies_z( coupler , state , state_tend , tracers , tracer_tend , space_op.recon , space_op.hydrostasis , dt );
        apply_tendencies             ( state , state_tend , tracers , tracer_tend , dt );
      } else {
        space_op.compute_tendencies_z( coupler , state , state_tend , tracers , tracer_tend , space_op.recon , space_op.hydrostasis , dt );
        apply_tendencies             ( state , state_tend , tracers , tracer_tend , dt );

        space_op.compute_tendencies_y( coupler , state , state_tend , tracers , tracer_tend , space_op.recon , dt );
        apply_tendencies             ( state , state_tend , tracers , tracer_tend , dt );

        space_op.compute_tendencies_x( coupler , state , state_tend , tracers , tracer_tend , space_op.recon , dt );
        apply_tendencies             ( state , state_tend , tracers , tracer_tend , dt );
      }
      dim_switch = ! dim_switch;

      #ifdef PAM_DEBUG
        std::vector<real> mass_final = compute_mass( coupler , state , tracers );
        for (int l=0; l < mass_final.size(); l++) {
          real mass_diff;
          if (mass_init[l] > 0) {
            mass_diff = abs(mass_final[l] - mass_init[l]) / abs(mass_init[l]);
          } else {
            mass_diff = mass_final[l];
          }
          real tol = 1.e-12;
          if (std::is_same<real,float>::value) {tol = 1.e-5;}
          if (mass_diff > tol) {
            std::cout << "Dycore mass change is too large. Abs Diff: " << abs(mass_final[l] - mass_init[l])
                      << ";   Rel Diff: " << mass_diff
                      << ";   Initial Mass: " << mass_init[l] << std::endl;
            // endrun("ERROR: mass not conserved by dycore");
          }
        }
        validate_array_positive(tracers);
        validate_array_inf_nan(state);
        validate_array_inf_nan(tracers);
      #endif

      if (coupler.get_num_dycore_functions() > 0) {
        awfl::tendencies_rho_theta::convert_dynamics_to_coupler_state( coupler , state , tracers );
        coupler.run_dycore_functions( dt );
        awfl::tendencies_rho_theta::convert_coupler_state_to_dynamics( coupler , state , tracers , space_op.hydrostasis );
      }

    }

    awfl::tendencies_rho_theta::convert_dynamics_to_coupler_state( coupler , state , tracers );
  }


  void finalize(pam::PamCoupler &coupler) { }


  const char * dycore_name() const { return "AWFL"; }

};
