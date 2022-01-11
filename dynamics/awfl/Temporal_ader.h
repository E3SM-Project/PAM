
#pragma once

#include "awfl_const.h"
#include "pam_coupler.h"

using pam::PamCoupler;

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = true;
int  constexpr nAder       = ngll;

template <class Spatial> class Temporal_operator {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  int  sponge_cells;
  real sponge_strength;

  real5d stateTend;
  real5d tracerTend;

  Spatial space_op;

  void init(int ny, int nx, int nens, real xlen, real ylen, int num_tracers, PamCoupler &coupler) {
    space_op.init(ny, nx, nens, xlen, ylen, num_tracers, coupler);

    std::string inFile = coupler.get_note( "standalone_input_file" );

    if (inFile != "") {
      YAML::Node config = YAML::LoadFile(inFile);
      sponge_cells    = config["sponge_cells"   ].as<int>();
      sponge_strength = config["sponge_strength"].as<real>();
    } else {
      sponge_cells    = 0;
      sponge_strength = 0;
    }

    stateTend  = space_op.createStateTendArr ();
    tracerTend = space_op.createTracerTendArr();
  }


  void convert_dynamics_to_coupler_state( DataManager &dm ) {
    space_op.convert_dynamics_to_coupler_state( dm );
  }


  void convert_coupler_state_to_dynamics( DataManager &dm ) {
    space_op.convert_coupler_state_to_dynamics( dm );
  }


  int add_tracer(DataManager &dm , std::string name , std::string desc , bool pos_def , bool adds_mass) {
    return space_op.add_tracer(dm , name , desc , pos_def , adds_mass);
  }


  void init_idealized_state_and_tracers( PamCoupler &coupler ) {
    space_op.init_idealized_state_and_tracers( coupler );
  }


  void output(PamCoupler &coupler, real etime) {
    space_op.output(coupler , etime);
  }


  real compute_time_step(PamCoupler &coupler, real cfl_in = -1) {
    return space_op.compute_time_step(coupler, cfl_in);
  }


  std::vector<real> compute_mass( PamCoupler &coupler , real5d &state , real5d &tracers ) {
    int nz = coupler.dm.get_dimension_size("z");
    int ny = coupler.dm.get_dimension_size("y");
    int nx = coupler.dm.get_dimension_size("x");
    int nens = coupler.dm.get_dimension_size("nens");

    int idR = Spatial::idR;
    int hs  = Spatial::hs;
    int num_tracers = space_op.num_tracers;
    YAKL_SCOPE( dz          , space_op.dz          );
    YAKL_SCOPE( hyDensCells , space_op.hyDensCells );

    std::vector<real> mass(num_tracers+1);
    real4d tmp("tmp",nz,ny,nx,nens);

    parallel_for( "Temporal_ader.h state mass" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tmp(k,j,i,iens) = (state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens)) * dz(k,iens);
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


  void timeStep( PamCoupler &coupler , real dtphys ) {
    YAKL_SCOPE( stateTend       , this->stateTend           );
    YAKL_SCOPE( tracerTend      , this->tracerTend          );
    YAKL_SCOPE( sponge_cells    , this->sponge_cells        );
    YAKL_SCOPE( sponge_strength , this->sponge_strength     );
    YAKL_SCOPE( hyDensCells     , this->space_op.hyDensCells);

    real dt = compute_time_step( coupler );

    real5d state   = space_op.createStateArr();
    real5d tracers = space_op.createTracerArr();

    #ifdef PAM_DEBUG
      memset(state   , 0._fp);
      memset(tracers , 0._fp);
    #endif

    space_op.convert_coupler_state_to_dynamics( coupler , state , tracers );

    int idR         = space_op.idR;
    int idU         = space_op.idU;
    int idV         = space_op.idV;
    int idW         = space_op.idW;
    int idT         = space_op.idT;
    int nx          = space_op.nx;
    int ny          = space_op.ny;
    int nz          = space_op.nz;
    int nens        = space_op.nens;
    int num_state   = space_op.num_state;
    int num_tracers = space_op.num_tracers;
    int hs          = space_op.hs;

    int n_iter = ceil( dtphys / dt );
    dt = dtphys / n_iter;

    for (int iter = 0; iter < n_iter; iter++) {

      #ifdef PAM_DEBUG
        validate_array_positive(tracers);
        validate_array_inf_nan(state);
        validate_array_inf_nan(tracers);
        std::vector<real> mass_init = compute_mass( coupler , state , tracers );
      #endif

      ScalarLiveOut<bool> neg_too_large(false);

      // Loop over different items in the spatial splitting
      for (int spl = 0 ; spl < space_op.numSplit() ; spl++) {
        real dtloc = dt;

        // Compute the tendencies for state and tracers
        space_op.computeTendencies( state , stateTend , tracers , tracerTend , dtloc , spl );

        parallel_for( "Temporal_ader.h apply tendencies" , SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          for (int l=0; l < num_state; l++) {
            state(l,hs+k,hs+j,hs+i,iens) += dtloc * stateTend(l,k,j,i,iens);
          }
          for (int l=0; l < num_tracers; l++) {
            tracers(l,hs+k,hs+j,hs+i,iens) += dtloc * tracerTend(l,k,j,i,iens);
            #ifdef PAM_DEBUG
              if (tracers(l,hs+k,hs+j,hs+i,iens) < -1.e-10) {
                neg_too_large = true;
              }
            #endif
            tracers(l,hs+k,hs+j,hs+i,iens) = max( 0._fp , tracers(l,hs+k,hs+j,hs+i,iens) );
          }
        });
      }

      #ifdef PAM_DEBUG
        if (neg_too_large.hostRead()) {
          std::cerr << "WARNING: Correcting a non-machine-precision negative tracer value" << std::endl;
          // endrun();
        }
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

      space_op.switch_directions();

      if (sponge_cells > 0) {
        real2d zint = coupler.dm.get<real,2>("vertical_interface_height");
        real2d zmid = coupler.dm.get<real,2>("vertical_midpoint_height");
        parallel_for( "Sponge" , SimpleBounds<4>(sponge_cells,ny,nx,nens) , YAKL_LAMBDA (int kk, int j, int i, int iens) {
          int k = nz-1-kk;
          real z1 = zint(nz-sponge_cells,iens);
          real z2 = zint(nz             ,iens);
          real znorm = (zmid(k,iens)-z1) / (z2 - z1);
          real mult = 1 - sponge_strength * ( cos(M_PI*znorm - M_PI) + 1 ) * 0.5_fp;
          real hydens = hyDensCells(k,iens);
          real dens_old = state(idR,hs+k,hs+j,hs+i,iens) + hydens;
          // state(idR,hs+k,hs+j,hs+i,iens) *= mult;
          // state(idU,hs+k,hs+j,hs+i,iens) *= mult;
          // state(idV,hs+k,hs+j,hs+i,iens) *= mult;
          state(idW,hs+k,hs+j,hs+i,iens) *= mult;
          state(idT,hs+k,hs+j,hs+i,iens) *= mult;
          // real dens_new = state(idR,hs+k,hs+j,hs+i,iens) + hydens;
          // for (int tr=0; tr < num_tracers; tr++) {
          //   real trac = tracers(tr,hs+k,hs+j,hs+i,iens);
          //   trac = trac / dens_old;
          //   trac *= mult;
          //   tracers(tr,hs+k,hs+j,hs+i,iens) = trac * dens_new;
          // }
        });
      }

    }

    space_op.convert_dynamics_to_coupler_state( coupler , state , tracers );
  }


  void finalize(PamCoupler &coupler) { }


  const char * dycore_name() const { return "AWFL"; }

};
