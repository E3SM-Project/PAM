
#pragma once

#include "const.h"
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
  
  void init(std::string inFile, int num_tracers, DataManager &dm) {
    space_op.init(inFile, num_tracers, dm);
    stateTend  = space_op.createStateTendArr ();
    tracerTend = space_op.createTracerTendArr();

    YAML::Node config = YAML::LoadFile(inFile);

    this->nens = config["nens"].as<int>();
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


  template <class MICRO>
  void timeStep( DataManager &dm , MICRO const &micro , real dt ) {
    YAKL_SCOPE( stateTend  , this->stateTend  );
    YAKL_SCOPE( tracerTend , this->tracerTend );

    real5d state   = dm.get<real,5>("dyn_state");
    real5d tracers = dm.get<real,5>("dyn_tracers");

    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < space_op.numSplit() ; spl++) {
      real dtloc = dt;

      // Compute the tendencies for state and tracers
      space_op.computeTendencies( state , stateTend , tracers , tracerTend , micro , dtloc , spl );

      int nx          = space_op.nx;
      int ny          = space_op.ny;
      int nz          = space_op.nz;
      int num_state   = space_op.num_state;
      int num_tracers = space_op.num_tracers;
      int hs          = space_op.hs;

      parallel_for( SimpleBounds<4>(nz,ny,nx,nens) , YAKL_LAMBDA (int k, int j, int i, int iens) {
        for (int l=0; l < num_state; l++) {
          state  (l,hs+k,hs+j,hs+i,iens) += dtloc * stateTend (l,k,j,i,iens);
        }
        for (int l=0; l < num_tracers; l++) {
          tracers(l,hs+k,hs+j,hs+i,iens) += dtloc * tracerTend(l,k,j,i,iens);
        }
      });
    }
  }


  void finalize(DataManager &dm) { }


  const char * getTemporalName() const { return "ADER-DT"; }

};

