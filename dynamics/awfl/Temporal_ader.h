
#pragma once

#include "const.h"

int  constexpr nTimeDerivs = 1;
bool constexpr timeAvg     = true;
int  constexpr nAder       = ngll;

/* REQUIRED:
int  static constexpr nTimeDerivs = [# tendency time derivatives needed by the time stepping scheme];
bool static constexpr timeAvg     = [whether the spatial operator should return a time-averaged tendency];

static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

void init(std::string inFile):
    - Process input file with filename "inFile"
    - Allocate and initialize internal stuff

void timeStep( StateArr &state , real dt ): 
    - Perform a single time step

const char * getTemporalName():
    - Return the name and info about this temporal operator
*/
template <class Spatial> class Temporal_ader {
public:
  static_assert(nTimeDerivs <= ngll , "ERROR: nTimeDerivs must be <= ngll.");

  typedef typename Spatial::StateTendArr  StateTendArr;
  typedef typename Spatial::StateArr      StateArr;
  typedef typename Spatial::TracerTendArr TracerTendArr;
  typedef typename Spatial::TracerArr     TracerArr;
  typedef typename Spatial::Location      Location;

  StateTendArr  stateTendArr;
  TracerTendArr tracerTendArr;

  Spatial spaceOp;
  
  void init(std::string inFile, int numTracers) {
    spaceOp.init(inFile, numTracers);
    stateTendArr  = spaceOp.createStateTendArr ();
    tracerTendArr = spaceOp.createTracerTendArr();
  }


  void timeStep( StateArr &stateArr , TracerArr &tracerArr , real dt ) {
    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < spaceOp.numSplit() ; spl++) {
      real dtloc = dt;
      spaceOp.computeTendencies( stateArr , stateTendArr , tracerArr , tracerTendArr , dtloc , spl );

      {
        auto &stateTendArr  = this->stateTendArr ;
        auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
          real &state      = spaceOp.getState     (stateArr     ,loc  ,spl);
          real &stateTend  = spaceOp.getStateTend (stateTendArr ,loc,0,spl);
          state  += dtloc * stateTend;
        };
        spaceOp.applyStateTendencies( applySingleTendency , spl );
      }

      {
        auto &tracerTendArr  = this->tracerTendArr ;
        auto applySingleTendency = YAKL_LAMBDA (Location const &loc) {
          real &tracer     = spaceOp.getTracer    (tracerArr    ,loc  ,spl);
          real &tracerTend = spaceOp.getTracerTend(tracerTendArr,loc,0,spl);
          tracer += dtloc * tracerTend;
        };
        spaceOp.applyTracerTendencies( applySingleTendency , spl );
      }
    }
  }


  void finalize(StateArr &state, TracerArr &tracers) {
    spaceOp.finalize( state , tracers );
  }


  const char * getTemporalName() const { return "ADER-DT"; }

};

