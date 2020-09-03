
#pragma once

#include "const.h"
#include "DataManager.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = true;
int  constexpr nAder       = ngll;

template <class Spatial> class Temporal_ader {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  typedef typename Spatial::Location      Location;

  real4d stateTend;
  real4d tracerTend;

  Spatial spaceOp;
  
  void init(std::string inFile, int num_tracers, DataManager &dm) {
    spaceOp.init(inFile, num_tracers, dm);
    stateTend  = spaceOp.createStateTendArr ();
    tracerTend = spaceOp.createTracerTendArr();
  }


  int add_tracer(DataManager &dm , std::string name , std::string desc , bool pos_def , bool adds_mass) {
    return spaceOp.add_tracer(dm , name , desc , pos_def , adds_mass);
  }


  template <class MICRO>
  void init_state( DataManager &dm , MICRO const &micro ) {
    spaceOp.init_state( dm , micro );
  }


  template <class F, class MICRO>
  void init_tracer_by_location(std::string name , F const &init_mass , DataManager &dm, MICRO const &micro) const {
    spaceOp.init_tracer_by_location(name , init_mass , dm, micro);
  }


  template <class MICRO>
  void adjust_state_for_moisture(DataManager &dm , MICRO const &micro) const {
    spaceOp.adjust_state_for_moisture( dm , micro );
  }


  template <class MICRO>
  real compute_time_step(real cfl, DataManager &dm, MICRO const &micro) {
    return spaceOp.compute_time_step(cfl, dm, micro);
  }


  template <class MICRO>
  void timeStep( DataManager &dm , MICRO const &micro , real dt ) {
    real4d state   = spaceOp.createStateArr();
    real4d tracers = spaceOp.createTracerArr();
    spaceOp.read_state_and_tracers( dm , state , tracers );

    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      real dtloc = dt;

      // Compute the tendencies for state and tracers
      spaceOp.computeTendencies( state , stateTend , tracers , tracerTend , micro , dtloc , spl );

      int nx          = spaceOp.nx;
      int ny          = spaceOp.ny;
      int nz          = spaceOp.nz;
      int num_state   = spaceOp.num_state;
      int num_tracers = spaceOp.num_tracers;
      int hs          = spaceOp.hs;

      parallel_for( Bounds<3>(nz,ny,nx) , YAKL_LAMBDA (int k, int j, int i) {
        for (int l=0; l < num_state; l++) {
          state(l,hs+k,hs+j,hs+i) += dtloc * stateTend(l,k,j,i);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,hs+i) += dtloc * tracerTend(l,k,j,i);
        }
      });
    }

    spaceOp.write_state_and_tracers( dm , state , tracers );
  }


  void finalize(DataManager &dm) { }


  const char * getTemporalName() const { return "ADER-DT"; }

};

