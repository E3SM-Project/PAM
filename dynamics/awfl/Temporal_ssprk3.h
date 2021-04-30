
#pragma once

#include "awfl_const.h"
#include "DataManager.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = false;
int  constexpr nAder       = 1;

template <class Spatial> class Temporal_operator {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  real5d stateTmp;
  real5d tracersTmp;

  real5d stateTend;
  real5d tracerTend;

  real5d stateTendAccum;
  real5d tracerTendAccum;

  Spatial space_op;
  
  void init(std::string inFile, int ny, int nx, int nens, real xlen, real ylen, int num_tracers, DataManager &dm) {
    space_op.init(inFile, ny, nx, nens, xlen, ylen, num_tracers, dm);
    stateTmp        = space_op.createStateArr     ();
    tracersTmp      = space_op.createTracerArr    ();

    stateTend       = space_op.createStateTendArr ();
    tracerTend      = space_op.createTracerTendArr();

    stateTendAccum  = space_op.createStateTendArr ();
    tracerTendAccum = space_op.createTracerTendArr();
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


  inline void zero_accum_arrays( real5d &stateTendAccum , real5d &tracerTendAccum ) {
    int num_state   = stateTendAccum .dimension[0];
    int nz          = stateTendAccum .dimension[1];
    int ny          = stateTendAccum .dimension[2];
    int nx          = stateTendAccum .dimension[3];
    int nens        = stateTendAccum .dimension[4];

    int num_tracers = tracerTendAccum.dimension[0];

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        stateTendAccum(l,k,j,i,iens) = 0;
      }
      for (int l=0; l < num_tracers; l++) {
        tracerTendAccum(l,k,j,i,iens) = 0;
      }
    });
  }


  inline void tendency_accum( real5d &stateTendAccum  , real5d const &stateTend ,
                              real5d &tracerTendAccum , real5d const &tracerTend ) {
    int num_state   = stateTendAccum .dimension[0];
    int nz          = stateTendAccum .dimension[1];
    int ny          = stateTendAccum .dimension[2];
    int nx          = stateTendAccum .dimension[3];
    int nens        = stateTendAccum .dimension[4];

    int num_tracers = tracerTendAccum.dimension[0];

    parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
      for (int l=0; l < num_state; l++) {
        stateTendAccum(l,k,j,i,iens) += stateTend(l,k,j,i,iens);
      }
      for (int l=0; l < num_tracers; l++) {
        tracerTendAccum(l,k,j,i,iens) += tracerTend(l,k,j,i,iens);
      }
    });
  }


  template <class MICRO>
  void timeStep( DataManager &dm , MICRO const &micro , real dtphys ) {
    YAKL_SCOPE( stateTend       , this->stateTend       );
    YAKL_SCOPE( stateTmp        , this->stateTmp        );
    YAKL_SCOPE( stateTendAccum  , this->stateTendAccum  );
    YAKL_SCOPE( tracerTend      , this->tracerTend      );
    YAKL_SCOPE( tracersTmp      , this->tracersTmp      );
    YAKL_SCOPE( tracerTendAccum , this->tracerTendAccum );

    real5d state   = dm.get<real,5>("dynamics_state");
    real5d tracers = dm.get<real,5>("dynamics_tracers");

    int nx          = space_op.nx;
    int ny          = space_op.ny;
    int nz          = space_op.nz;
    int nens        = space_op.nens;
    int num_state   = space_op.num_state;
    int num_tracers = space_op.num_tracers;
    int hs          = space_op.hs;

    real dt = compute_time_step( 0.8 , dm , micro );

    real loctime = 0.;
    while (loctime < dtphys) {
      if (loctime + dt > dtphys) { dt = dtphys - loctime; }

      real dtloc = dt;

      /////////////////////////////////////
      // Stage 1
      /////////////////////////////////////
      zero_accum_arrays( stateTendAccum , tracerTendAccum );
      for (int spl = 0 ; spl < space_op.numSplit() ; spl++) {
        space_op.computeTendencies( state , stateTend , tracers , tracerTend , micro , dtloc , spl );
        tendency_accum( stateTendAccum , stateTend , tracerTendAccum , tracerTend );
      }
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          stateTmp  (l,hs+k,hs+j,hs+i,iens) = state  (l,hs+k,hs+j,hs+i,iens) + dtloc * stateTendAccum (l,k,j,i,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracersTmp(l,hs+k,hs+j,hs+i,iens) = tracers(l,hs+k,hs+j,hs+i,iens) + dtloc * tracerTendAccum(l,k,j,i,iens);
        }
      });

      /////////////////////////////////////
      // Stage 2
      /////////////////////////////////////
      zero_accum_arrays( stateTendAccum , tracerTendAccum );
      for (int spl = 0 ; spl < space_op.numSplit() ; spl++) {
        space_op.computeTendencies( stateTmp , stateTend , tracersTmp , tracerTend , micro , dtloc , spl );
        tendency_accum( stateTendAccum , stateTend , tracerTendAccum , tracerTend );
      }
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          stateTmp(l,hs+k,hs+j,hs+i,iens) = 0.75_fp * state   (l,hs+k,hs+j,hs+i,iens) + 
                                            0.25_fp * stateTmp(l,hs+k,hs+j,hs+i,iens) +
                                            0.25_fp * dtloc * stateTendAccum(l,k,j,i,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracersTmp(l,hs+k,hs+j,hs+i,iens) = 0.75_fp * tracers   (l,hs+k,hs+j,hs+i,iens) + 
                                              0.25_fp * tracersTmp(l,hs+k,hs+j,hs+i,iens) +
                                              0.25_fp * dtloc * tracerTendAccum(l,k,j,i,iens);
        }
      });

      /////////////////////////////////////
      // Stage 3
      /////////////////////////////////////
      zero_accum_arrays( stateTendAccum , tracerTendAccum );
      for (int spl = 0 ; spl < space_op.numSplit() ; spl++) {
        space_op.computeTendencies( stateTmp , stateTend , tracersTmp , tracerTend , micro , dtloc , spl );
        tendency_accum( stateTendAccum , stateTend , tracerTendAccum , tracerTend );
      }
      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          state(l,hs+k,hs+j,hs+i,iens) = (1._fp/3._fp) * state   (l,hs+k,hs+j,hs+i,iens) + 
                                         (2._fp/3._fp) * stateTmp(l,hs+k,hs+j,hs+i,iens) +
                                         (2._fp/3._fp) * dtloc * stateTendAccum(l,k,j,i,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,hs+i,iens) = (1._fp/3._fp) * tracers   (l,hs+k,hs+j,hs+i,iens) + 
                                           (2._fp/3._fp) * tracersTmp(l,hs+k,hs+j,hs+i,iens) +
                                           (2._fp/3._fp) * dtloc * tracerTendAccum(l,k,j,i,iens);
        }
      });

      loctime += dt;
    }
  }


  void finalize(DataManager &dm) { }


  const char * dycore_name() const { return "AWFL"; }

};

