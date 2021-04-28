
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = true;
int  constexpr nAder       = ngll;

template <class Spatial> class Temporal_operator {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  int nens;

  real5d stateTend;
  real5d tracerTend;

  Spatial space_op;
  
  void init(std::string inFile, int ny, int nx, int nens, real xlen, real ylen, int num_tracers, DataManager &dm) {
    space_op.init(inFile, ny, nx, nens, xlen, ylen, num_tracers, dm);
    stateTend  = space_op.createStateTendArr ();
    tracerTend = space_op.createTracerTendArr();
  }


  template <class MICRO>
  void convert_dynamics_to_coupler_state( DataManager &dm , MICRO &micro ) {
    space_op.convert_dynamics_to_coupler_state( dm , micro );
  }


  template <class MICRO>
  void convert_coupler_state_to_dynamics( DataManager &dm , MICRO &micro ) {
    space_op.convert_coupler_state_to_dynamics( dm , micro );
  }


  int add_tracer(DataManager &dm , std::string name , std::string desc , bool pos_def , bool adds_mass) {
    return space_op.add_tracer(dm , name , desc , pos_def , adds_mass);
  }


  template <class MICRO>
  void init_state_and_tracers( DataManager &dm , MICRO const &micro ) {
    space_op.init_state_and_tracers( dm , micro );
  }


  template <class F, class MICRO>
  void init_tracer_by_location(std::string name , F const &init_mass , DataManager &dm, MICRO const &micro) const {
    space_op.init_tracer_by_location(name , init_mass , dm, micro);
  }


  template <class MICRO>
  void output(DataManager &dm, MICRO const &micro, real etime) const {
    space_op.output(dm , micro , etime);
  }


  template <class MICRO>
  real compute_time_step(real cfl, DataManager &dm, MICRO const &micro) {
    return space_op.compute_time_step(cfl, dm, micro);
  }


  std::vector<real> compute_mass( DataManager &dm ) {
    real5d state   = dm.get<real,5>("dynamics_state");
    real5d tracers = dm.get<real,5>("dynamics_tracers");
    int nz = dm.get_dimension_size("z");
    int ny = dm.get_dimension_size("y");
    int nx = dm.get_dimension_size("x");
    int nens = dm.get_dimension_size("nens");

    int idR = Spatial::idR;
    int hs  = Spatial::hs;
    int num_tracers = space_op.num_tracers;
    YAKL_SCOPE( dz          , space_op.dz          );
    YAKL_SCOPE( hyDensCells , space_op.hyDensCells );

    std::vector<real> mass(num_tracers+1);
    real4d tmp("tmp",nz,ny,nx,nens);

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      tmp(k,j,i,iens) = (state(idR,hs+k,hs+j,hs+i,iens) + hyDensCells(k,iens)) * dz(k,iens);
    });
    mass[0] = yakl::intrinsics::sum(tmp);

    for (int l=0; l < num_tracers; l++) {
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        tmp(k,j,i,iens) = tracers(l,hs+k,hs+j,hs+i,iens) * dz(k,iens);
      });
      mass[l+1] = yakl::intrinsics::sum(tmp);
    }
    return mass;
  }


  template <class MICRO>
  void timeStep( DataManager &dm , MICRO const &micro , real dtphys ) {
    YAKL_SCOPE( stateTend  , this->stateTend  );
    YAKL_SCOPE( tracerTend , this->tracerTend );

    real5d state   = dm.get<real,5>("dynamics_state");
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    real dt = compute_time_step( 0.8 , dm , micro );

    real loctime = 0.;
    while (loctime < dtphys) {
      if (loctime + dt > dtphys) { dt = dtphys - loctime; }

      #ifdef PAM_DEBUG
        validate_array_positive(tracers);
        validate_array_inf_nan(state);
        validate_array_inf_nan(tracers);
        std::vector<real> mass_init = compute_mass( dm );
      #endif

      ScalarLiveOut<bool> neg_too_large(false);

      // Loop over different items in the spatial splitting
      for (int spl = 0 ; spl < space_op.numSplit() ; spl++) {
        real dtloc = dt;

        // Compute the tendencies for state and tracers
        space_op.computeTendencies( state , stateTend , tracers , tracerTend , micro , dtloc , spl );

        int nx          = space_op.nx;
        int ny          = space_op.ny;
        int nz          = space_op.nz;
        int nens        = space_op.nens;
        int num_state   = space_op.num_state;
        int num_tracers = space_op.num_tracers;
        int hs          = space_op.hs;

        parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
          for (int l=0; l < num_state; l++) {
            state(l,hs+k,hs+j,hs+i,iens) += dtloc * stateTend(l,k,j,i,iens);
          }
          for (int l=0; l < num_tracers; l++) {
            tracers(l,hs+k,hs+j,hs+i,iens) += dtloc * tracerTend(l,k,j,i,iens);
            if (tracers(l,hs+k,hs+j,hs+i,iens) < -1.e-10) {
              neg_too_large = true;
            }
            tracers(l,hs+k,hs+j,hs+i,iens) = max( 0._fp , tracers(l,hs+k,hs+j,hs+i,iens) );
          }
        });
      }

      #ifdef PAM_DEBUG
        if (neg_too_large.hostRead()) {
          std::cerr << "WARNING: Correcting a non-machine-precision negative tracer value" << std::endl;
          // endrun();
        }
        std::vector<real> mass_final = compute_mass( dm );
        for (int l=0; l < mass_final.size(); l++) {
          real mass_diff;
          if (mass_init[l] > 0) {
            mass_diff = abs(mass_final[l] - mass_init[l]) / abs(mass_init[l]);
          } else {
            mass_diff = mass_final[l];
          }
          if (mass_diff > 1.e-12) {
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

      loctime += dt;
    }
  }


  void finalize(DataManager &dm) { }


  const char * getTemporalName() const { return "ADER-DT"; }

};

